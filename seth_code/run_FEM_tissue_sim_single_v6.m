
% Units
%
% voltage in mV
% current in uA
% conductance in mS
% resistance in k-ohm
% capacitance in uF    (uF*mV/ms = uA)
% time in ms
% length in um
% concentration in mM
%
% New in v3: added option for single cell simulations to determine
% steady-state for tissue simulation initial conditions (using specified
% bcl and current scaling factors)

s = 'r.-'; clear p;
model = 'LR1';      % guinea pig ventricle
%model = 'Court98';  % human atrial
%model = 'TNNP04';  % human ventricle
%model = 'ORd11';    % human ventricle
%model = 'ORd11_Clancy'; % w/Na mutant channels from Clancy papers
% model = 'Bond04';  % mouse ventricle
% model = 'Bond04_Clancy';  % mouse ventricle w/Na mutant channels
%model = 'LRd07_Clancy'; % guinea pig w/Na mutant channels

IDdist = 'chan'; % 'area'; %'specified'; % distribution within the ID patches

FEM_file1 = 'FEMDATA_V_100Parts_IP16nm_P16nm_Map1_Gjy_Wvy_MJ_0_peri_0_chan.mat';
FEM_file2 = 'FEMDATA_V_100Parts_IP60nm_P60nm_Map1_Gjy_Wvy_MJ_0_peri_0_chan.mat';

clamp_Na = 0;
clamp_K = 0;

Mdisc = 100;

scale_INa = 1;
scale_IKr = 1;
mf1 = 0;
mf2 = 1;

single_cell_ic_flag = 0; % option to run single cell simulations to obtain new initial conditions, otherwise default ICs
% for future: better version would be to run simulations for varying bcl
% and scaleI -> save/load ICs
cable_cv_est_flag = 0;
transmural_flag = 0;

Ncell = 20;     % number of cells
% cell geometry parameters
L = 100;        % cell length, um
r = 5.5;         % cell radius, um

% estimate propagation time down cable, based on a set conduction velocity value, 
% used to determine twin (see below)
% for future: could also use a short monodomain simulation to estimate CV
cv_est = 20*1e3/100;  % um/ms (D = 0.1 ~ 20 cm/s, use 10 cm/s)
Lcable = L*Ncell; % um
tprop_est = Lcable/cv_est;  % ms
% time parameters
bcl = 20;      % basic cycle length ms
nbeats = 1;      % number of beats
T = bcl*nbeats;

% model-specific time steps
switch model
    case {'Bond04','Bond04_Clancy'}
        dt_scale = .01;
    case 'ORd11_Clancy'
        dt_scale = 0.2; 
    case 'LRd07_Clancy'
        dt_scale = 0.1;
    otherwise
        dt_scale = 1;
end

% different numerical integration and sampling methods
if 0
    dtparams.method = 'fixed';  % 'fixed', 'adapt', or 'adapt_samp'
    dtparams.dt1 = dt_scale*1e-2;      % ms
    dtparams.dt1_samp = dtparams.dt1*1; 
    dtparams.dt2 = dtparams.dt1*50;   
    dtparams.dt2_samp = dtparams.dt2*1;  
    dtparams.dtS1 = dtparams.dt1*0.2; % time step for cleft concentration integration
    dtparams.dtS2 = dtparams.dt1*1;
    dtparams.twin = max(25, round(tprop_est)+10); % defines window using dt1
elseif 0
   dtparams.method = 'adapt';
   dtparams.dt_min = dt_scale*1e-2;  % min time step used during upstroke (< twin), ms
   dtparams.dt_max = dtparams.dt_min*50;  % max time step, ms
   dtparams.dtS_min = dtparams.dt_min*0.2; % min time step for cleft concentration integration, used while dt = dt_min, should be <= dt_min
   dtparams.dtS_max = dtparams.dt_min*2;  % max time step for cleft concentration, must be <= dt_max
   dtparams.dV_max = 0.5;  % max voltage change, mV
   dtparams.dS_max = 0.01;  % max cleft concentration change, mM
   dtparams.twin = max(25, round(tprop_est)+10); % defines window using dt_min
else
   dtparams.method = 'adapt_samp';
   dtparams.dt_min = dt_scale*1e-2;  % min time step used during upstroke (< twin), ms
   dtparams.dt_min_samp = dtparams.dt_min/dt_scale; % sampling time
   dtparams.dt_max = dtparams.dt_min*10;  % max time step, ms
   dtparams.dt_max_samp = dtparams.dt_max*10/dt_scale;
   dtparams.dtS_min = dtparams.dt_min*1; % min time step for cleft concentration integration, used while dt = dt_min, should be <= dt_min
   dtparams.dtS_max = dtparams.dt_min*2;  % max time step for cleft concentration, must be <= dt_max
   dtparams.dV_max = 0.5;  % max voltage change, mV
   dtparams.dS_max = 0.01;  % max cleft concentration change, mM
%    dtparams.twin = max(25, round(tprop_est)+10); % defines window using dt_min
   dtparams.twin = max(5, round(tprop_est)+2); % defines window using dt_min

end
Vthresh = -60;  % activation/repol threshold, mV
single_beat_flag = 0; % 0 to stop simulation after propagation of a single beat
trange = [max(0,(nbeats-2)*bcl-10) nbeats*bcl];  % time range to output values
%trange = [0 nbeats*bcl];  % time range to output values


% model-independent parameters

% load FEM data
cd ../mesh_data; 
load(FEM_file1); FEM_data1 = FEM_data;
load(FEM_file2); FEM_data2 = FEM_data;
cd ../seth_code; 

p_ext = 150*10;  % extracellular resistivity, k-ohm*um
f_disc = 1; f_bulk = 1; % cleft conductance scaling factors
fVol_cleft = 1; % cleft volume scaling factor

Gc_array = zeros(Mdisc, Mdisc, Ncell);
Gb_mat = zeros(Ncell, Mdisc);
IDarea_vec = zeros(Mdisc, 2*Ncell);
chan_area_norm_mat = zeros(Mdisc, 2*Ncell);
Vol_cleft_vec = zeros(Mdisc, Ncell-1);

% define chan_area_norm term (if not in FEM_data file)
FEM_data1.chan_area_norm = FEM_data1.partition_surface/sum(FEM_data1.partition_surface);
FEM_data2.chan_area_norm = FEM_data2.partition_surface/sum(FEM_data2.partition_surface);

% note junction i is between pre-junctional membrane on cell i and
% post-junctional membrane on cell i+1
ind_junc = [8 14];  % junction for FEM_data2 file properties (all others use FEM_data1)
% set baseline for all cells/junctions
for i = 1:Ncell 
    Gc_array(:,:,i) = f_disc*(FEM_data1.cleft_adjacency_matrix)/p_ext;  % mS, Mdisc x Mdisc
    Gb_mat(i,:) = f_bulk*FEM_data1.bulk_adjacency_matrix'/p_ext;  % mS, 1 x Mdisc
    IDarea_vec(:,2*i-1) = FEM_data1.partition_surface;  % ID membrane patch surface area, um^2
    IDarea_vec(:,2*i) = FEM_data1.partition_surface;  % ID membrane patch surface area, um^2
    chan_area_norm_mat(:,2*i-1) = FEM_data1.chan_area_norm;
    chan_area_norm_mat(:,2*i) = FEM_data1.chan_area_norm;
    if i < Ncell
        Vol_cleft_vec(:,i) = fVol_cleft*FEM_data1.partition_volume; % cleft volume, um^3
    end  
end
for i = 1:Ncell-1 % loop over all junctions
    if ismember(i,ind_junc)
        Gc_array(:,:,i) = f_disc*(FEM_data2.cleft_adjacency_matrix)/p_ext;  % mS, Mdisc x Mdisc
        Gb_mat(i,:) = f_bulk*FEM_data2.bulk_adjacency_matrix'/p_ext;  % mS, 1 x Mdisc
        IDarea_vec(:,2*i) = FEM_data2.partition_surface;  % ID membrane patch surface area, um^2
        IDarea_vec(:,2*i+1) = FEM_data2.partition_surface;  % ID membrane patch surface area, um^2
        chan_area_norm_mat(:,2*i) = FEM_data2.chan_area_norm;
        chan_area_norm_mat(:,2*i+1) = FEM_data2.chan_area_norm;
        Vol_cleft_vec(:,i) = fVol_cleft*FEM_data2.partition_volume; % cleft volume, um^3       
    end
end

% gap junctional coupling
D = 1*ones(Ncell-1,1);      % diffusion coefficient, cm^2/s, used to determine ggap
% D(round(Ncell/2)) = 1;  % vary one junction in the cable middle

% cell geometry
Aax = 2*pi*r*L; % patch surface area, um^2
% add some variation
Ad = sum(IDarea_vec);
Atot = zeros(1, Ncell);
for i = 1:Ncell
    Atot(i) = Ad(2*i-1) + Ad(2*i) + Aax; % total surface area, um^2, 1 x Ncell vector
end
% bulk extracellular concentrations
Ko = 5.4;                  % mM
Nao = 140;                 % mM
Cao = 1.8;                 % mM

% ID localization (vs axial patch)
locUniform = (Ad(1:2:end-1)+Ad(2:2:end))/2./Atot;  % value for uniform distribution (proportional to area)
locINa = 0.9;  % overall localization at ID
locIK1 = 0.5;
switch model
    case 'LR1'
        clear locINaK;
    case 'ORd11'
        locINaK = 0.1;
end

% ID distribution (must sum to 1)
if strcmp(IDdist,'specified') % must have size Mdisc x 1
    IDdistINa = [.4; .3; .2; .1];  % distrubtion of channels in ID nanodomains
    IDdistIK1 = ones(Mdisc,1)/Mdisc;
    IDdistINaK = ones(Mdisc,1)/Mdisc;
end

Npatches = 2*Mdisc*(Ncell-1) + Ncell; % number of membrane patches

% storage parameters
storeparams.int_store_all = 100;  % ms, interval for storing all state variables
% gating variable indices to store (model specific), [ ] = store all variables, NaN = store no variables
storeparams.indG_store = nan;


% model specific parameters
switch model
    case 'LR1'
        [p, x0] = InitialConstants_LR91(Atot');
        
        p.Isi_max = constructCell2Patches(Ncell, Mdisc, p.Isi_max);
        p.IKp_max = constructCell2Patches(Ncell, Mdisc, p.IKp_max);
        p.Ib_max = constructCell2Patches(Ncell, Mdisc, p.Ib_max);
        p.INa_max = constructCell2Patches(Ncell, Mdisc, p.INa_max);
        p.IK1_max = constructCell2Patches(Ncell, Mdisc, p.IK1_max);
        p.IK_max = constructCell2Patches(Ncell, Mdisc, p.IK_max);
        % order is determined by code in fun_name
        % INa, Isi, IK, IK1, IKp, Ib
        Ncurrents = 6;
        scaleI = ones(1,Ncurrents);
        
        p.iina = 1; p.iisi = 2; p.iik = 3;
        p.iik1 = 4; p.iikp = 5; p.iib = 6;
        loc_vec = zeros(1, Ncurrents); % ID localization vec
        loc_vec(p.iina) = locINa;
        loc_vec(p.iik1) = locIK1;
        
        loc_vec(p.iik) = 0;
        loc_vec(p.iikp) = 0;
        
        scaleI(p.iik1) = 1;
        
        % ionic model-specific parameters
        ionic_fun_name = 'fun_LR1';
        
        % initial conditions
        p.Nstate = 8;  % number of state variables, including Vm, per patch
        p.mLR1 = 1; % flag for modified LR1 model with Ca2+ speedup
        
        % stimulus parameters
        p.stim_dur = 1;   % ms
        p.stim_amp = 80e-8*mean(Atot);    % uA
        p.istim = ((1:5)-1)*(2*Mdisc+1)+1;  % indices of axial membrane patch
        p.indstim = ismember((1:Npatches)',p.istim);
        
    case 'ORd11'
        
        % order is determined by code in fun_name
        % INa INaL Ito ICaL IKr IKs IK1 INaCa_i INaCa_ss INaK  IKb INab ICab IpCa
        Ncurrents = 14;
        scaleI = ones(1, Ncurrents); % ionic current scaling factors
        %          scaleI(2) = 5; scaleI(5) = .15; % generates EADs for bcl = 1000, homogeneous endo, D=1 cable
        
        p.iina = 1; p.iinal = 2; p.iito = 3;
        p.iical = 4; p.iikr = 5; p.iiks = 6;
        p.iik1 = 7; p.iinaca_i = 8; p.iinaca_ss = 9;
        p.iinak = 10; p.iikb = 11; p.iinab = 12;
        p.iicab = 13; p.iipca = 14;
        
        scaleI(p.iik1) = 3;
        scaleI(p.iina) = 3;
        
        loc_vec = zeros(1, Ncurrents); % ID localization vec
        loc_vec(p.iina) = locINa;
        loc_vec(p.iik1) = locIK1;
        loc_vec(p.iinak) = locINaK;
        
        loc_vec(p.iito) = mean(locUniform);
        loc_vec(p.iikr) = mean(locUniform);
        loc_vec(p.iiks) = mean(locUniform);
        
        % additional scaling factors
        p.fSERCA = 1; p.fRyR = 1; p.ftauhL = 1;
        p.fCaMKa = 1; p. fIleak = 1; p.fJrel = 1;
        
        % ionic model-specific parameters
        ionic_fun_name = 'fun_ORd11';
        
        % initial conditions
        %initial conditions for state variables
        x0 = Initial_ORd11;
%         x0 = 1e2*[-0.879989146999539, 0.070944539841272, 0.070945359018008, 1.448320027716102, 1.448319784601420, 0.000000851590737, 0.000000840540955, 0.016028979649317, 0.015570572412294, 0.000073463295678, 0.006980036544213, 0.006979869883167, 0.006978935065780, 0.004548302689923, 0.006978307799113, 0.000001883686661, 0.005011493681533, 0.002694910834916, 0.000010012992091, 0.009995539443664, 0.005900573961157, 0.000005101890248, 0.009995539516859, 0.006425686271562, 0.000000000023424, 0.009999999908676, 0.009093809723890, 0.009999999908675, 0.009998130840372, 0.009999751803949, 0.000026423594695, 0.009999999908620, 0.009999999908648, 0.000000080785953, 0.004517981715199, 0.002727141574241, 0.000001929116438, 0.009967605243606, 0.000000002422481, 0.000000003026518, 0.000122276951936];
        %X0 is the vector for initial sconditions for state variables
        
        p.Nstate = 41;  % number of state variables, including Vm, per patch
        
        % stimulus parameters
        p.stim_dur = 1;   % ms
        p.stim_amp = 50;    % uA/uF
        stim_cells = 1:5;  % cell indices for stimulus
        p.istim = (stim_cells-1)*(2*Mdisc+1)+1;  % indices of axial membrane patches
        p.indstim = ismember((1:Npatches)',p.istim);
        
        % cell type
        p.celltype = zeros(Npatches,1); %endo = 0, epi = 1, M = 2
        if transmural_flag
            % transmural wedge: endo -> M -> epi
            % Cells 1:iEndoM -> endo
            % Cells iEndoM+1:iMEpi -> M
            % Cells iMEpi+1:end -> epi
            iEndoM = 60; iMEpi = 105;
            p.celltype(1+iEndoM+(2*iEndoM-1)*Mdisc:iMEpi+(2*iMEpi-1)*Mdisc) = 2;
            p.celltype(1+iMEpi+(2*iMEpi-1)*Mdisc:end) = 1;
        end
        
    case 'ORd11_Clancy'
        % order is determined by code in fun_name
        % 1   2    3   4    5   6   7   8       9        10    11  12   13   14
        % INa INaL Ito ICaL IKr IKs IK1 INaCa_i INaCa_ss INaK  IKb INab ICab IpCa
        Ncurrents = 14;
        scaleI = ones(1, Ncurrents); % ionic current scaling factors
        %         scaleI(1) = 1.5; scaleI(2) = 10; scaleI(4) = 2; scaleI(5) = .5;
    
        
        p.iina = 1; p.iinal = 2; p.iito = 3;
        p.iical = 4; p.iikr = 5; p.iiks = 6;
        p.iik1 = 7; p.iinaca_i = 8; p.iinaca_ss = 9;
        p.iinak = 10; p.iikb = 11; p.iinab = 12;
        p.iicab = 13; p.iipca = 14;
        
        scaleI(p.iina) = scale_INa; 
        scaleI(p.iinal) = 0; 
        scaleI(p.iical) = 1;
        scaleI(p.iikr) = scale_IKr;
        
        loc_vec = zeros(1, Ncurrents); % ID localization vec
        loc_vec(p.iina) = locINa;
%         loc_vec(p.iinal) = locINaL;
        loc_vec(p.iik1) = locIK1;
%         loc_vec(p.iinak) = locINaK;
        
        % additional scaling factors
        p.fSERCA = 1; p.fRyR = 1; p.ftauhL = 1;
        p.fCaMKa = 1; p. fIleak = 1; p.fJrel = 1;
        
        % ionic model-specific parameters
        ionic_fun_name = 'fun_ORd11_Clancy';
        
        % mutant-specific parameters
        p.MCflag = 1;  % 1 = Monte Carlo Na model, 0 = HH gating Na model
        % heterozygous channels - two populations
        % mutant_flag: 0 or 1 (0 for WT)
        % mutant_type: 'Y1795C', 'I1768V', 'dKPQ', '1795insD'
        p.mutant_flag1 = mf1; p.mutant_type1 = 'Y1795C';
        p.mutant_flag2 = mf2; p.mutant_type2 = 'Y1795C';
        
        x0 = Initial_ORd11_Clancy(p.mutant_type1, p.mutant_flag1, p.mutant_type2, p.mutant_flag2);
        %X0 is the vector for initial sconditions for state variables
        
        p.Nstate = 67;  % number of state variables, including Vm, per patch
        
        % stimulus parameters
        p.stim_dur = 1;   % ms
        p.stim_amp = 50;    % uA/uF
        p.istim = ((1:5)-1)*(2*Mdisc+1)+1;  % indices of axial membrane patches
        p.indstim = ismember((1:Npatches)',p.istim);
        
        % cell type
        p.celltype = zeros(Npatches,1); %endo = 0, epi = 1, M = 2
        if 0
            % transmural wedge: endo -> M -> epi
            % Cells 1-iEndoM -> endo
            % Cells iEndoM+1:iMEpi -> M
            % Cells iMEpi+1:end -> epi
            iEndoM = 15; iMEpi = 25;
            p.celltype(1+iEndoM+(2*iEndoM-1)*Mdisc:iMEpi+(2*iMEpi-1)*Mdisc) = 2;
            p.celltype(1+iMEpi+(2*iMEpi-1)*Mdisc:end) = 1;
        end
    case 'Court98'
        % order is determined by code in fun_name
        % INa IK1 Ito IKur IKr IKs IBNa IBK IBCa INaK ICaP INaCa ICaL
        Ncurrents = 13;
        scaleI = ones(1, Ncurrents); % ionic current scaling factors
        
        p.iina = 1; p.iik1 = 2; p.iito = 3;
        p.iikur = 4; p.iikr = 5; p.iiks = 6;
        p.iibna = 7; p.iibk = 8; p.iibca = 9;
        p.iinak = 10; p.iicap = 11; p.iinaca = 12;
        p.iical = 13;
        
        loc_vec = zeros(1, Ncurrents); % ID localization vec
        loc_vec(p.iina) = locINa;
        %         loc_vec(p.iik1) = locIK1;
        %         loc_vec(p.iinak) = locINaK;
        %
        %         loc_vec(p.iito) = locUniform;
        %         loc_vec(p.iikr) = locUniform;
        %         loc_vec(p.iiks) = locUniform;
        %
        
        % ionic model-specific parameters
        ionic_fun_name = 'fun_courtemachne98';
        
        % initial conditions
        %initial conditions for state variables
        x0 = Initial_Court98;
        %X0 is the vector for initial sconditions for state variables
        
        p.Nstate = 21;  % number of state variables, including Vm, per patch
        
        % stimulus parameters
        p.stim_dur = 1;   % ms
        p.stim_amp = -60;    % uA/uF
        stim_cells = 1:5;  % cell indices for stimulus
        p.istim = (stim_cells-1)*(2*Mdisc+1)+1;  % indices of axial membrane patches
        p.indstim = ismember((1:Npatches)',p.istim);
        
    case 'TNNP04'
        % order is determined by code in fun_name
        % INa IK1 Ito IKr IKs IBNa IpK IBCa INaK ICaP INaCa ICaL
        Ncurrents = 12;
        scaleI = ones(1, Ncurrents); % ionic current scaling factors
        
        p.iina = 1; p.iik1 = 2; p.iito = 3;
        p.iikr = 4; p.iiks = 5; p.iibna = 6;
        p.iipk = 7; p.iibca = 8; p.iinak = 9;
        p.iipca = 10; p.iinaca = 11; p.iical = 12;
        
        loc_vec = zeros(1, Ncurrents); % ID localization vec
        loc_vec(p.iina) = locINa;
        %         loc_vec(p.iik1) = locIK1;
        %         loc_vec(p.iinak) = locINaK;
        %
        %         loc_vec(p.iito) = locUniform;
        %         loc_vec(p.iikr) = locUniform;
        %         loc_vec(p.iiks) = locUniform;
        %
        
        % ionic model-specific parameters
        ionic_fun_name = 'fun_tnnp04';
        
        % initial conditions
        %initial conditions for state variables
        x0 = Initial_TNNP04;
        %X0 is the vector for initial sconditions for state variables
        
        p.Nstate = 17;  % number of state variables, including Vm, per patch
        
        % stimulus parameters
        p.stim_dur = 1;   % ms
        p.stim_amp = -52;    % uA/uF
        stim_cells = 1:5;  % cell indices for stimulus
        p.istim = (stim_cells-1)*(2*Mdisc+1)+1;  % indices of axial membrane patches
        p.indstim = ismember((1:Npatches)',p.istim);
        
        % cell type
        p.celltype = zeros(Npatches,1); %endo = 0, epi = 1, M = 2
        
    case 'Bond04'
        % order is determined by code in fun_name
        % ICaL IpCa INaCa ICaB INa INab INaK IKto_f IKto_s IK1 IKs IKur IKss IKr IClCa
        Ncurrents = 15;
        scaleI = ones(1, Ncurrents); % ionic current scaling factors
        
        p.iical = 1; p.iipca = 2; p.iinaca = 3;
        p.iicab = 4; p.iina = 5; p.iinab = 6;
        p.iinak = 7; p.iikto_f = 8; p.iikto_s = 9;
        p.iik1 = 10; p.iiks = 11; p.iikur = 12;
        p.iikss = 13; p.iikr = 14; p.iiclca = 15;
        
        loc_vec = zeros(1, Ncurrents); % ID localization vec
        loc_vec(p.iina) = locINa;
        %         loc_vec(p.iik1) = locIK1;
        %         loc_vec(p.iinak) = locINaK;
        %
        %         loc_vec(p.iito) = locUniform;
        %         loc_vec(p.iikr) = locUniform;
        %         loc_vec(p.iiks) = locUniform;
        %
        
        % ionic model-specific parameters
        ionic_fun_name = 'fun_bondarenko04';
        
        % initial conditions
        %initial conditions for state variables
        x0 = Initial_Bondarenko04;
        %X0 is the vector for initial sconditions for state variables
        
        p.Nstate = 41;  % number of state variables, including Vm, per patch
        
        % stimulus parameters
        p.stim_dur = 1;   % ms
        p.stim_amp = -40;    % uA/uF
        stim_cells = 1:5;  % cell indices for stimulus
        p.istim = (stim_cells-1)*(2*Mdisc+1)+1;  % indices of axial membrane patches
        p.indstim = ismember((1:Npatches)',p.istim);
        
         % cell type
        p.celltype = zeros(Npatches,1); % apical = 0, septum = 1
        
    case 'Bond04_Clancy'
            % order is determined by code in fun_name
        % ICaL IpCa INaCa ICaB INa INab INaK IKto_f IKto_s IK1 IKs IKur
        % IKss IKr IClCa
        Ncurrents = 15;
        scaleI = ones(1, Ncurrents); % ionic current scaling factors
        
        p.iical = 1; p.iipca = 2; p.iinaca = 3;
        p.iicab = 4; p.iina = 5; p.iinab = 6;
        p.iinak = 7; p.iikto_f = 8; p.iikto_s = 9;
        p.iik1 = 10; p.iiks = 11; p.iikur = 12;
        p.iikss = 13; p.iikr = 14; p.iiclca = 15;
       
        scaleI(p.iikur) = .75;

        % ionic model-specific parameters
        ionic_fun_name = 'fun_bondarenko04_Clancy';
        
        % mutant-specific parameters
        p.MCflag = 0;  % Markov chain Na model, 0 = original model, 1 = Clancy models
        % heterozygous channels - two populations
        % mutant_flag: 0 or 1 (0 for WT)
        % mutant_type: 'Y1795C', 'I1768V', 'dKPQ', '1795insD'
        p.mutant_flag1 = 0; p.mutant_type1 = 'dKPQ';
        p.mutant_flag2 = 1; p.mutant_type2 = 'dKPQ';
        
        if p.mutant_flag1 || p.mutant_flag2
            p.MCflag = 1;
        end
        
        x0 = Initial_Bondarenko04_Clancy(p.mutant_type1, p.mutant_flag1, p.mutant_type2, p.mutant_flag2);
        %X0 is the vector for initial sconditions for state variables
        
        %X0 is the vector for initial sconditions for state variables
        
        p.Nstate = 67;  % number of state variables, including Vm, per patch
        
         % stimulus parameters
        p.stim_dur = 1;   % ms
        p.stim_amp = -40;    % uA/uF
        stim_cells = 1:5;  % cell indices for stimulus
        p.istim = (stim_cells-1)*(2*Mdisc+1)+1;  % indices of axial membrane patches
        p.indstim = ismember((1:Npatches)',p.istim);
        
         % cell type
        p.celltype = zeros(Npatches,1); % apical = 0, septum = 1
        
    case 'LRd07_Clancy'
        
        Ncurrents = 10;
        scaleI = ones(1, Ncurrents); % ionic current scaling factors
              
        % order is determined by code in fun_name
        % INa, INab, IKs, IKr, IK1, IKp, ICaL, ICaT, INaCa, INaK
        p.iina = 1; p.iinab = 2; p.iiks = 3;
        p.iikr = 4; p.iik1 = 5; p.iikp = 6; 
        p.iical = 7; p.iicat = 8; p.iinaca = 9; p.iinak = 10;
        
        % initial cleft ion concentrations
        loc_vec = zeros(1, Ncurrents); % ID localization vec
        loc_vec(p.iina) = locINa;
        loc_vec(p.iik1) = locIK1;
%         loc_vec(p.iinak) = locINaK;
        
        % ionic model-specific parameters
        ionic_fun_name = 'fun_LRd07_Clancy';
        
        % mutant-specific parameters
        p.MCflag = 1;  % 1 = Monte Carlo Na model, 0 = HH gating Na model
        % mutant_flag: 0 or 1 (0 for WT)
        % mutant_type: 'Y1795C', 'I1768V', 'dKPQ', '1795insD'
        p.mutant_flag = 1; p.mutant_type = 'dKPQ';
        p.Atot = Atot;     
        [p2,x0] = InitialConstants_LRd07_Clancy(r*1e-4, L*1e-4, p.mutant_type, p.mutant_flag);
        %X0 is the vector for initial conditions for state variables
        
        mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
        p = mergestructs(p,p2);
        
        p.Nstate = 33;  % number of state variables, including Vm, per patch
        
        % stimulus parameters
        p.stim_dur = 1;   % ms
        p.stim_amp = 1.5*-50;    % uA/uF
        p.istim = ((1:5)-1)*(2*Mdisc+1)+1;  % indices of axial membrane patches
        p.indstim = ismember((1:Npatches)',p.istim);
        
      
end


switch IDdist
    case 'equal'
        % equal distribution in the ID nanodomains
        loc_mat = repmat(loc_vec, Mdisc, 1)/Mdisc;  % Mdisc x Ncurrents, ID localization matrix
    case 'random'
        % random distribution (uniform)
        tmp = rand(Mdisc,1); tmp = tmp/sum(tmp);
        loc_mat = tmp*loc_vec;
    case 'prop'
        % proportional to cleft width
        tmp = w_vec; tmp = tmp/sum(tmp);
        loc_mat = tmp*loc_vec;
    case 'inverseprop'
        % inversely proportional to cleft
        tmp = 1./w_vec; tmp = tmp/sum(tmp);
        loc_mat = tmp*loc_vec;
    case 'specified'
        loc_mat = zeros(Mdisc, length(loc_vec)); % Mdisc x Ncurrents, ID localization matrix
        loc_mat(:,p.iina) = IDdistINa*locINa;
        loc_mat(:,p.iik1) = IDdistIK1*locIK1;
        if exist('locINaK','var')
            loc_mat(:,p.iinak) = IDdistINaK*locINaK;
        end
    case 'area'
        loc_mat = zeros(2*Mdisc, Ncurrents, Ncell);
        for i = 1:Ncell
            tmp = IDarea_vec(:,2*i-1); tmp = tmp/sum(tmp);
            loc_mat(1:Mdisc,:,i) = tmp*loc_vec/2;
             
            tmp = IDarea_vec(:,2*i); tmp = tmp/sum(tmp);
            loc_mat(1+Mdisc:2*Mdisc,:,i) = tmp*loc_vec/2;
        end
    case 'chan'
        loc_mat = zeros(2*Mdisc, Ncurrents, Ncell);
        % for channels with area-based distribution
        for i = 1:Ncell
            tmp = IDarea_vec(:,2*i-1); tmp = tmp/sum(tmp);
            loc_mat(1:Mdisc,:,i) = tmp*loc_vec/2;
            
            tmp = IDarea_vec(:,2*i); tmp = tmp/sum(tmp);
            loc_mat(1+Mdisc:2*Mdisc,:,i) = tmp*loc_vec/2;
        end
        ind = [p.iina p.iik1];  % channel-based distribution
        for i = 1:Ncell
            tmp = chan_area_norm_mat(:,2*i-1); tmp = tmp/sum(tmp);
            loc_mat(1:Mdisc,ind,i) = tmp*loc_vec(ind)/2;
            
            tmp = chan_area_norm_mat(:,2*i); tmp = tmp/sum(tmp);
            loc_mat(1+Mdisc:2*Mdisc,ind,i) = tmp*loc_vec(ind)/2;
        end
    case 'region_spec'
        
        loc_ID = zeros(Mdisc, Ncurrents); % Mdisc x Ncurrents, localization matrix within ID
        
        locInterpNa = 0.7;  % fraction of ID Na channels in interplicate
        locInterpK1 = nan;  % fraction of ID K1 channels in interplicate
        % convention: nan = uniform by area
        
        % identify plicate / interplicate regions
        zpos = FEM_data.partition_centers(:,3);
        plicate_thresh = [0.2 0.9]; % z-position defining plicate / interplicate
        ind_plicate = find(zpos < plicate_thresh(1) | zpos > plicate_thresh(2));
        ind_interplicate = setdiff(1:Mdisc, ind_plicate);
        
        % fraction of area of plicate / interplicate regions (for reference)
        totArea = sum(IDarea_vec);
        plicateArea = sum(IDarea_vec(ind_plicate));
        interplicateArea = sum(IDarea_vec(ind_interplicate));
        frac_plicateArea = plicateArea/totArea;
        frac_interplicateArea = interplicateArea/totArea;
        
        locInterplicate = zeros(Ncurrents,1);  % fraction of channels in interplicate (out of all ID)
        locInterplicate(p.iina) = locInterpNa;  % defined 
        locInterplicate(p.iik1) = locInterpK1;
        locInterplicate(isnan(locInterplicate)) = frac_interplicateArea;  % proportional to area
        locPlicate = 1 - locInterplicate;
        
        tmp = IDarea_vec(ind_interplicate); tmp = tmp/sum(tmp);
        loc_ID(ind_interplicate,:) = tmp*locInterplicate';
        tmp = IDarea_vec(ind_plicate); tmp = tmp/sum(tmp);
        loc_ID(ind_plicate,:) = tmp*locPlicate';
        
        loc_mat = loc_ID.*(ones(Mdisc,1)*loc_vec);  % Mdisc x Ncurrents, ID localization matrix
end


p.L = L; p.r = r;  % um
% extracellular concentrations
p.K_o = Ko;                  % mM
p.Na_o = Nao;                 % mM
p.Ca_o = Cao;                 % mM
p.A_o = Ko + Nao + 2*Cao;
% cleft ionic concentration parameters
p.cleft_conc = [p.Na_o; p.K_o; p.Ca_o; p.A_o];
clamp_flag = [clamp_Na; clamp_K; 0; 0];  % Na, K, Ca (1 to clamp cleft conc, 0 for dynamic cleft concentration)
zvec = [1; 1; 2; -1];  % charge valence

if single_cell_ic_flag % option to run single cell simulations to obtain new initial conditions, fixed extracellular ionic concentrations
    nbeat_single = 10;
    % use fixed two time step method for single cell sims
    dt1 = 1e-2; dt2 = 1e-1; twin = 25; 
    psingle = p; psingle.dt = dt1; psingle.bcl = bcl;
    psingle.Npatches = 1; psingle.f_I = scaleI'; 
    psingle.indstim = 1; 
    psingle.celltype = 0;
    psingle.ind_tau_ip = 1; psingle.ind_tau_im = 1; psingle.tau_Nai_mat = 1; psingle.tau_Cai_mat = 1; psingle.tau_Ki_mat = 1;
    Cm = 1*1e-8;      % membrane capacitance, uF/um^2
    psingle.Ctot = (2*sum(IDarea_vec) + 2*pi*r*L)*Cm;   % total cell capacitance, uF
    Vm_single(1) = x0(1); Gsingle = x0(2:end)';
    ionic_fun = str2func(['@(t,x,p,S) ',ionic_fun_name,'(t,x,p,S)']);
    
    ti = 0;  
    while ti < nbeat_single*bcl
        if mod(ti,bcl)<twin
            dt = dt1; 
        else
            dt = dt2; 
        end
        psingle.dt = dt;
        [Gsingle, Iion, ~] = ionic_fun(ti,[Vm_single;Gsingle],psingle, p.cleft_conc);
        Vm_single = Vm_single - dt*Iion/psingle.Ctot;
        
        ti = round(ti + dt,4);  
    end
    x0 = [Vm_single; Gsingle]; 
end

if cable_cv_est_flag
    pcable = p; dtparams_cable = dtparams;
    switch dtparams.method
        case 'fixed'
            dtparams_cable.dt2 = dtparams_cable.dt1; 
        case {'adapt','adapt_samp'}
            dtparams_cable.dt_max = dtparams_cable.dt_min;
    end
    pcable.phii_0 = x0(1);
    pcable.indstim = ismember((1:Ncell)',1:5);
    pcable.stim_amp = pcable.stim_amp/2;
    pcable.ind_tau_ip = 1; pcable.ind_tau_im = 1; pcable.tau_Nai_mat = 1; pcable.tau_Cai_mat = 1; pcable.tau_Ki_mat = 1;
    g0_vec = zeros(Ncell*(pcable.Nstate-1),1);
    for i = 1:pcable.Nstate-1
        g0_vec(1 + Ncell*(i-1):Ncell*i) = x0(i+1);
    end
    pcable.g0_vec = g0_vec;
    pcable.celltype = zeros(Ncell,1); %endo = 0, epi = 1, M = 2
    if transmural_flag
        % transmural wedge: endo -> M -> epi
        % Cells 1:iEndoM -> endo
        % Cells iEndoM+1:iMEpi -> M
        % Cells iMEpi+1:end -> epi
        iEndoM = 60; iMEpi = 105;
        pcable.celltype(1+iEndoM:iMEpi) = 2;
        pcable.celltype(1+iMEpi:end) = 1;
    end
    
    [~, ~, ~, ~, ~, cv_est] =  generic_cable(ionic_fun_name, ...
        pcable, Ncell, D, bcl, scaleI, dtparams_cable, T, Vthresh, 1, trange); % CV in cm/s  
    tprop_est = 2 * Lcable/(cv_est*1e3/100);  % ms
    % Note: monodomain cable sims typically produce CV values greater than
    % the ephaptic model.  Include scaling factor to account for this
    dtparams.twin = max(25, round(tprop_est)+15); % defines window using dt1

end

p.phii_0 = x0(1);
g0_vec = zeros(Npatches*(p.Nstate-1),1);
for i = 1:p.Nstate-1
    g0_vec(1 + Npatches*(i-1):Npatches*i) = x0(i+1);
end
p.g0_vec = g0_vec;

[ts, phi_mat, Gmat, Smat, Iall_mat, tup, trepol, pout, cv, store_out] =  generic_IDnano_epc_v6(ionic_fun_name, ...
    p, Ncell, Mdisc, Gb_mat, Gc_array, D, bcl, loc_mat, IDarea_vec, scaleI, Vol_cleft_vec, clamp_flag, ...
    zvec, dtparams, T, Vthresh, single_beat_flag, trange, storeparams);

% indices
Nphi = 3*Ncell-2 + Mdisc*(Ncell-1);  % number of voltage state variables
Ncleft_comp = (Ncell-1)*Mdisc;  % number of cleft compartments
Npatches = 2*Mdisc*(Ncell-1) + Ncell; % number of membrane patches

ind_Vm = ((1:Ncell)-1)*(2*Mdisc+1)+1;  % indices of axial membrane patch
ind_phim = ((1:Ncell)-1)*(Mdisc+3)+1; % indices of axial potentials
ind_phiR = ((1:Ncell-1)-1)*(Mdisc+3)+2; % indices of right potentail
ind_phiL =  ((2:Ncell)-1)*(Mdisc+3);  % indices of left potentials
ind_phic = setdiff(1:Nphi, [ind_phim ind_phiR ind_phiL]);  % indices of cleft potentials
ind_VR = []; ind_VL = [];
for i = 1:Ncell-1, for j = 1:Mdisc, ind_VR = [ind_VR, (i-1)*(2*Mdisc+1)+j+1]; end, end
for i = 2:Ncell, for j = 1:Mdisc, ind_VL = [ind_VL, (i-2)*(2*Mdisc+1)+Mdisc+j+1]; end, end

icell = 2;  % indices for a specific cell
i_phii_cell = icell;
i_VL_cell = (icell-2)*Mdisc+1:(icell-2)*Mdisc+Mdisc;
i_VR_cell = (icell-1)*Mdisc+1:(icell-1)*Mdisc+Mdisc;

Vm_mat = constructVm_mat(phi_mat, Npatches, Mdisc, Ncell);

phi_i = phi_mat(ind_phim,:);
phi_cleft = phi_mat(ind_phic,:);
phi_L = phi_mat(ind_phiL,:);
phi_R = phi_mat(ind_phiR,:);

% note: Vm = phi_i = Vm_mat(ind_Vm,:)
VL = Vm_mat(ind_VL,:);
VR = Vm_mat(ind_VR,:);
Na_cleft = Smat(1:Ncleft_comp,:);
K_cleft = Smat(Ncleft_comp+1:2*Ncleft_comp,:);
Ca_cleft = Smat(2*Ncleft_comp+1:3*Ncleft_comp,:);

if 1
    figure(1); dd = 5;  % note this only plots a sub-set of positions
  
    st = 'krbgmck'; count = 0;
    for q = 1:dd:Ncell
        count = count + 1;
        
        subplot(2,2,1); plot(ts, phi_i(q,:),s); hold on;
        subplot(2,2,2); plot(ts, phi_cleft((q-1)*Mdisc+1:q*Mdisc,:),s); hold on;
        
        subplot(2,2,3); plot(ts, Na_cleft((q-1)*Mdisc+1:q*Mdisc,:),s); hold on;
        subplot(2,2,4); plot(ts, K_cleft((q-1)*Mdisc+1:q*Mdisc,:),s); hold on;
    end
    figure(2);
    plot(ts, phi_i(i_phii_cell,:),s); hold on;
%     plot(ts, VL(i_VL_cell,:),'r'); hold on;
%     plot(ts, VL(i_VR_cell,:),'b'); hold on;
    
    
        
    figure(3);  % plot ionic currents, uA (model specific)
    Ncurrents = length(scaleI);
    INa = Iall_mat(pout.iina:Ncurrents:end,:); 
    IK1 = Iall_mat(pout.iik1:Ncurrents:end,:);
    subplot(2,2,1);  plot(ts, INa); hold on;
    subplot(2,2,2);  plot(ts, IK1); hold on;
    
    % current densities (uA/uF)
    Cax = pout.Cax;
    CL = pout.CL(2:Ncell,:); CLvec = reshape(CL',(Ncell-1)*Mdisc,1);
    CR = pout.CR(1:Ncell-1,:); CRvec = reshape(CR',(Ncell-1)*Mdisc,1);
        
    Call_vec = nan(Npatches,1);
    Call_vec(ind_Vm) = Cax;
    Call_vec(ind_VL) = CLvec;
    Call_vec(ind_VR) = CRvec;
    
    subplot(2,2,3); plot(ts, INa./(Call_vec*ones(1,length(ts))));
    subplot(2,2,4); plot(ts, IK1./(Call_vec*ones(1,length(ts))));
    
%     icell = 10;
%     figure(3);
%     xlim = 1500+[2 3]; ylim = [min(phi_cleft(:)) max(phi_cleft(:))];
%     n = floor(sqrt(Mdisc));
%     for i = 1:Mdisc
%         subplot(n, n, i);
%         plot(ts, phi_cleft(icell*Mdisc+i,:),'k','linewidth',2);
%         set(gca,'xlim',xlim,'ylim',ylim); axis off;
%     end
%     
%     ind = find(ts>= (nbeats-1)*bcl+2 & ts<=(nbeats-1)*bcl+3.5);
%     t = ts(ind); 
%     tscale = 3;
%     Vscale = 300;
%     
%     generate_3d_mesh_plot3(FEM_data.partition_centers(:,1), FEM_data.partition_centers(:,2),...
%         FEM_data.partition_centers(:,3),FEM_data.cleft_adjacency_matrix, [.2 .8], 1);
%     
%     for idisc = 1:Mdisc
%         x = FEM_data.partition_centers(idisc,1);
%         y = FEM_data.partition_centers(idisc,2);
%         z = FEM_data.partition_centers(idisc,3);
%         V = phi_cleft(icell*Mdisc+idisc,ind);
%         
%         generate_3d_mesh_time_plot(x, y, z, t, V, tscale, Vscale,[.2 .8]); grid on;
%     end
%     set(gca,'fontsize',28);
%     
%     icell = 15; tval = 2505.65;
%     ind = find(ts>= (nbeats-1)*bcl & ts<=(nbeats-1)*bcl+10);
%     
%     subplot(4,3,1);
%     plot(ts, phi_i(icell,:),'k','linewidth',2);
%     hold on; plot(tval*[1 1],[-100 50],'k--','linewidth',2);
%     plot(2600+[0 100],[-80 -80],'k','linewidth',8);
%     set(gca,'xlim',[2500 2860],'fontsize',28,'box','off','tickdir','out','linewidth',2,'xtick',[],'xcolor','w');
% 
%     subplot(4,3,2);
%     plot(ts(ind), phi_i(icell,ind),'k','linewidth',2);
%     hold on; plot(tval*[1 1],[-100 50],'k--','linewidth',2);
%     plot(2505.9+[0 1],[-80 -80],'k','linewidth',8);
%     set(gca,'xlim',[2504 2507],'fontsize',28,'box','off','tickdir','out','linewidth',2,'xtick',[],'xcolor','w');
% 
%     subplot(4,3,3);
%     hold on; plot(tval*[1 1],[-60 5],'k--','linewidth',2);
%     plot(2505.9+[0 1],[-50 -50],'k','linewidth',8);
%     set(gca,'xlim',[2504 2507],'fontsize',28,'box','off','tickdir','out','linewidth',2,'xtick',[],'xcolor','w');
%     
%     for zzz = 1:2
%         switch zzz
%             case 1
%                 sp = 3:2:7;
%             case 2
%                 sp = 4:2:8;
%         end
%        
%         % ind = find(ts>= (nbeats-1)*bcl+2 & ts<=(nbeats-1)*bcl+20);
%         cvals = jet(1000);
%         crange = [-60 5];
%         %  crange = [90 140];
%         for idisc = 1:Mdisc
%             x = FEM_data.partition_centers(idisc,1);
%             y = FEM_data.partition_centers(idisc,2);
%             z = FEM_data.partition_centers(idisc,3);
%             V = phi_cleft(icell*Mdisc+idisc,ind);
%             s = cvals(round(interp1(linspace(crange(1),crange(2),1000),1:1000, min(V))),:);
%             
%             subplot(4,3,3);
%             plot(ts(ind), V,'color',[s .2]);hold on;
%             
%             subplot(4,2,sp);
%             plot3(x, y, z, '.','markersize',50,'color',s); hold on; grid on;
%             
%             
%         end
%         generate_3d_mesh_plot3(FEM_data.partition_centers(:,1), FEM_data.partition_centers(:,2),...
%             FEM_data.partition_centers(:,3),FEM_data.cleft_adjacency_matrix, [.2 .8], 0);
%         set(gca,'fontsize',28,'tickdir','in');
%         switch zzz
%             case 1
%                 view(-7, 89); set(gca,'xticklabel',[],'yticklabel',[],'zticklabel',[]); 
%                 hold on; plot3([8 10],[.5 .5],[0 0],'k','linewidth',8); grid on;
%             case 2
%                 colorbar; caxis(crange); colormap(jet);
%                 view(-20, 10); set(gca,'xticklabel',[],'yticklabel',[],'zticklabel',[]); 
%                 hold on; plot3([0 0],[8 10.1],[0 0],'k','linewidth',8);
%                 plot3([0 0],[10 10],[0 .2],'k','linewidth',8); grid on;
%               %  plot3([0 2],[10 10],[0 0],'k','linewidth',8);
%         end
%     end
%     
%       figure(4);
%      xlim = 1500+[0 30]; ylim = [min(Na_cleft(:)) max(Na_cleft(:))];
%     n = floor(sqrt(Mdisc));
%     for i = 1:Mdisc
%         subplot(n, n, i);
%         plot(ts, Na_cleft(icell*Mdisc+i,:),'k','linewidth',2);
%         set(gca,'xlim',xlim,'ylim',ylim); axis off;
%     end
%     
%     %
    if 0
        figure(1);
        v = VideoWriter('test.avi');
        open(v);
        ind = find(ts>= (nbeats-1)*bcl & ts<=(nbeats-1)*bcl+10);
        crange = [90 140];
        for i = 1:10:length(ind)
            V = Na_cleft(icell*Mdisc+1:(icell+1)*Mdisc,ind(i));
            plot_3d_mesh_voltage(FEM_data.partition_centers(:,1), FEM_data.partition_centers(:,2),...
                FEM_data.partition_centers(:,3), V, FEM_data.cleft_adjacency_matrix, crange);
            set(gca,'fontsize',28);
            title(['t = ',num2str(ts(ind(i)))]); grid on;
            hold off; view(165, -21);  colorbar; colormap(jet); caxis(crange);
            
            pause(.01); frame = getframe(gcf); writeVideo(v,frame);
        end
        close(v);
    end
    
    % analysis of ionic currents
    if 0
        % capacitance, uF
        Cax = pout.Cax;
        CL = pout.CL(2:Ncell,:); CLvec = reshape(CL',(Ncell-1)*Mdisc,1);
        CR = pout.CR(1:Ncell-1,:); CRvec = reshape(CR',(Ncell-1)*Mdisc,1);
        
        switch model
            case 'LR1'
                
                % INa
                ENa_ax =p.RTF*log(Nao./p.Na_i);       % Nernst potential of Na, mV
                ENa_cleft =p.RTF*log(Na_cleft./p.Na_i); % Nernst potential of Na, mV
                m = Gmat(1:Npatches,:);
                H = Gmat(1+Npatches:2*Npatches,:);
                J = Gmat(1+2*Npatches:3*Npatches,:);
                
                INa_ax = (pout.f_I(p.iina,ind_Vm)'*ones(1,length(ts))).*p.INa_max.*m(ind_Vm,:).^3.*H(ind_Vm,:).*J(ind_Vm,:).*(phi_i-ENa_ax);
                INa_L = (pout.f_I(p.iina,ind_VL)'*ones(1,length(ts))).*p.INa_max.*m(ind_VL,:).^3.*H(ind_VL,:).*J(ind_VL,:).*(VL-ENa_cleft);
                INa_R = (pout.f_I(p.iina,ind_VR)'*ones(1,length(ts))).*p.INa_max.*m(ind_VR,:).^3.*H(ind_VR,:).*J(ind_VR,:).*(VR-ENa_cleft);
                
                % IK1
                EK1_ax= p.RTF*log(Ko./p.K_i);
                ak1_ax = 1.02./(1+exp(0.2385*(phi_i-EK1_ax-59.215)));
                bk1_ax = (0.49124*exp(0.08032*(phi_i-EK1_ax+5.476))+exp(0.06175*(phi_i-EK1_ax-594.31)))./...
                    (1+exp(-0.5143*(phi_i-EK1_ax+4.753)));
                sq_ax=sqrt(Ko/5.4);
                gK1_ax =p.IK1_max*sq_ax.*ak1_ax./(ak1_ax+bk1_ax);
                IK1_ax = (pout.f_I(p.iik1,ind_Vm)'*ones(1,length(ts))).*gK1_ax.*(phi_i-EK1_ax);
                
                EK1_cleft= p.RTF*log(K_cleft./p.K_i);
                sq_cleft=sqrt(K_cleft/5.4);
                ak1_L = 1.02./(1+exp(0.2385*(VL-EK1_cleft-59.215)));
                bk1_L = (0.49124*exp(0.08032*(VL-EK1_cleft+5.476))+exp(0.06175*(VL-EK1_cleft-594.31)))./...
                    (1+exp(-0.5143*(VL-EK1_cleft+4.753)));
                gK1_L =p.IK1_max*sq_cleft.*ak1_L./(ak1_L+bk1_L);
                IK1_L = (pout.f_I(p.iik1,ind_VL)'*ones(1,length(ts))).*gK1_L.*(VL-EK1_cleft);
                
                ak1_R = 1.02./(1+exp(0.2385*(VR-EK1_cleft-59.215)));
                bk1_R = (0.49124*exp(0.08032*(VR-EK1_cleft+5.476))+exp(0.06175*(VR-EK1_cleft-594.31)))./...
                    (1+exp(-0.5143*(VR-EK1_cleft+4.753)));
                gK1_R =p.IK1_max*sq_cleft.*ak1_R./(ak1_R+bk1_R);
                IK1_R = (pout.f_I(p.iik1,ind_VR)'*ones(1,length(ts))).*gK1_R.*(VR-EK1_cleft);
                
            case 'ORd11'
                
                %physical constants
                R=8314.0; T=310.0; F=96485.0;  % C/mol
                Rcg = 2; % uF
                
                % INa
                nai = Gmat(1:Npatches,:);
                cass = Gmat(5*Npatches+1:6*Npatches,:);
                CaMKt = Gmat(39*Npatches+1:40*Npatches,:);
                m = Gmat(8*Npatches+1:9*Npatches,:);
                hf = Gmat(9*Npatches+1:10*Npatches,:);
                hs = Gmat(10*Npatches+1:11*Npatches,:);
                j = Gmat(11*Npatches+1:12*Npatches,:);
                hsp = Gmat(12*Npatches+1:13*Npatches,:);
                jp = Gmat(13*Npatches+1:14*Npatches,:);
                
                Ahf=0.99;Ahs=1.0-Ahf;
                h=Ahf*hf+Ahs*hs;
                hp=Ahf*hf+Ahs*hsp;
                GNa=75;
                %CaMK constants
                KmCaMK=0.15; CaMKo=0.05; KmCaM=0.0015;
                %update CaMK
                CaMKb=CaMKo*(1.0-CaMKt)./(1.0+KmCaM./cass);
                CaMKa=p.fCaMKa*(CaMKb+CaMKt);
                fINap=(1.0./(1.0+KmCaMK./CaMKa));
                RTF = (R*T/F);
                ENa_ax =RTF*log(Nao./nai(ind_Vm,:));       % Nernst potential of Na, mV
                ENa_L =RTF*log(Na_cleft./nai(ind_VL,:)); % Nernst potential of Na, mV
                ENa_R =RTF*log(Na_cleft./nai(ind_VR,:)); % Nernst potential of Na, mV
                
                % uA
                INa_ax = Rcg*pout.Ctot*(pout.f_I(p.iina,ind_Vm)'*ones(1,length(ts))).*GNa.*(phi_i-ENa_ax).*m(ind_Vm,:).^3.*((1-fINap(ind_Vm,:)).*h(ind_Vm,:).*j(ind_Vm,:)+fINap(ind_Vm,:).*hp(ind_Vm,:).*jp(ind_Vm,:));
                INa_L = Rcg*pout.Ctot*(pout.f_I(p.iina,ind_VL)'*ones(1,length(ts))).*GNa.*(VL-ENa_L).*m(ind_VL,:).^3.*((1-fINap(ind_VL,:)).*h(ind_VL,:).*j(ind_VL,:)+fINap(ind_VL,:).*hp(ind_VL,:).*jp(ind_VL,:));
                INa_R = Rcg*pout.Ctot*(pout.f_I(p.iina,ind_VR)'*ones(1,length(ts))).*GNa.*(VR-ENa_R).*m(ind_VR,:).^3.*((1-fINap(ind_VR,:)).*h(ind_VR,:).*j(ind_VR,:)+fINap(ind_VR,:).*hp(ind_VR,:).*jp(ind_VR,:));
                
                % IK1
                ki = Gmat(2*Npatches+1:3*Npatches,:);
                xk1 = Gmat(36*Npatches+1:37*Npatches,:);
                
                EK_ax =RTF*log(Ko./ki(ind_Vm,:));       % Nernst potential of K, mV
                EK_L =RTF*log(K_cleft./ki(ind_VL,:));
                EK_R =RTF*log(K_cleft./ki(ind_VR,:));
                
                rk1_ax = 1.0./(1.0+exp((phi_i+105.8-2.6*Ko)./9.493));
                rk1_L = 1.0./(1.0+exp((VL+105.8-2.6*K_cleft)./9.493));
                rk1_R = 1.0./(1.0+exp((VR+105.8-2.6*K_cleft)./9.493));
                
                GK1=0.1908*ones(Npatches,1); GK1(p.celltype==1) = 0.1908*1.2; GK1(p.celltype==2) = 0.1908*1.3;
                
                % uA
                IK1_ax = Rcg*pout.Ctot*((GK1(ind_Vm)'.*pout.f_I(p.iik1,ind_Vm))'*ones(1,length(ts))).*sqrt(Ko).*rk1_ax.*xk1(ind_Vm,:).*(phi_i-EK_ax);
                IK1_L = Rcg*pout.Ctot*((GK1(ind_VL)'.*pout.f_I(p.iik1,ind_VL))'*ones(1,length(ts))).*sqrt(K_cleft).*rk1_L.*xk1(ind_VL,:).*(VL-EK_L);
                IK1_R = Rcg*pout.Ctot*((GK1(ind_VR)'.*pout.f_I(p.iik1,ind_VR))'*ones(1,length(ts))).*sqrt(K_cleft).*rk1_R.*xk1(ind_VR,:).*(VR-EK_R);
                
                % INaK
                Nao_all = nan(Npatches,length(ts));
                Ko_all = nan(Npatches,length(ts));
                Nao_all(ind_Vm,:) = Nao;
                Ko_all(ind_Vm,:) = Ko;
                for ii = 1:Ncell-1
                    Nao_all(2+(ii-1)*(2*Mdisc+1):2+(ii-1)*(2*Mdisc+1)+Mdisc-1,:) = Na_cleft(1+(ii-1)*Mdisc:ii*Mdisc,:);
                    Nao_all(2+(ii-1)*(2*Mdisc+1)+Mdisc:ii*(2*Mdisc+1),:) = Na_cleft(1+(ii-1)*Mdisc:ii*Mdisc,:);
                    Ko_all(2+(ii-1)*(2*Mdisc+1):2+(ii-1)*(2*Mdisc+1)+Mdisc-1,:) = K_cleft(1+(ii-1)*Mdisc:ii*Mdisc,:);
                    Ko_all(2+(ii-1)*(2*Mdisc+1)+Mdisc:ii*(2*Mdisc+1),:) = K_cleft(1+(ii-1)*Mdisc:ii*Mdisc,:);
                end
                
                k1p=949.5; k1m=182.4; k2p=687.2; k2m=39.4; k3p=1899.0; k3m=79300.0; k4p=639.0;
                k4m=40.0; Knai0=9.073; Knao0=27.78; delta=-0.1550;
                Knai=Knai0*exp((delta*Vm_mat*F)/(3.0*R*T));
                Knao=Knao0*exp(((1.0-delta)*Vm_mat*F)/(3.0*R*T));
                Kki=0.5; Kko=0.3582; MgADP=0.05; MgATP=9.8; Kmgatp=1.698e-7; H=1.0e-7;
                eP=4.2; Khp=1.698e-7; Knap=224.0; Kxkur=292.0;
                P=eP./(1.0+H/Khp+nai./Knap+ki./Kxkur);
                a1=(k1p*(nai./Knai).^3.0)./((1.0+nai./Knai).^3.0+(1.0+ki./Kki).^2.0-1.0);
                b1=k1m*MgADP; a2=k2p; b2=(k2m*(Nao_all./Knao).^3.0)./((1.0+Nao_all./Knao).^3.0+(1.0+Ko_all./Kko).^2.0-1.0);
                a3=(k3p*(Ko_all/Kko).^2.0)./((1.0+Nao_all./Knao).^3.0+(1.0+Ko_all./Kko).^2.0-1.0);
                b3=(k3m*P.*H)./(1.0+MgATP./Kmgatp); a4=(k4p*MgATP/Kmgatp)./(1.0+MgATP./Kmgatp);
                b4=(k4m.*(ki./Kki).^2.0)./((1.0+nai./Knai).^3.0+(1.0+ki./Kki).^2.0-1.0);
                x1=a4.*a1.*a2+b2.*b4.*b3+a2.*b4.*b3+b3.*a1.*a2;
                x2=b2.*b1.*b4+a1.*a2.*a3+a3.*b1.*b4+a2.*a3.*b4;
                x3=a2.*a3.*a4+b3.*b2.*b1+b2.*b1.*a4+a3.*a4.*b1;
                x4=b4.*b3.*b2+a3.*a4.*a1+b2.*a4.*a1+b3.*b2.*a1;
                E1=x1./(x1+x2+x3+x4); E2=x2./(x1+x2+x3+x4);
                E3=x3./(x1+x2+x3+x4); E4=x4./(x1+x2+x3+x4);
                zk=1.0; zna=1.0; JnakNa=3.0*(E1.*a3-E2.*b3); JnakK=2.0*(E4.*b1-E3.*a1);
                Pnak=30*ones(Npatches,1); Pnak(p.celltype==1) = 30*0.9;
                Pnak(p.celltype==2) = 30*0.7;
                
                % uA
                INaK_ax = Rcg*pout.Ctot*((Pnak(ind_Vm)'.*pout.f_I(p.iinak,ind_Vm))'*ones(1,length(ts))).*(zna.*JnakNa(ind_Vm,:)+zk.*JnakK(ind_Vm,:));
                INaK_L = Rcg*pout.Ctot*((Pnak(ind_VL)'.*pout.f_I(p.iinak,ind_VL))'*ones(1,length(ts))).*(zna.*JnakNa(ind_VL,:)+zk.*JnakK(ind_VL,:));
                INaK_R = Rcg*pout.Ctot*((Pnak(ind_VR)'.*pout.f_I(p.iinak,ind_VR))'*ones(1,length(ts))).*(zna.*JnakNa(ind_VR,:)+zk.*JnakK(ind_VR,:));
                
        end
        
        % calculate densities, uA/uF
        INa_ax_density = INa_ax/Cax;
        INa_L_density = INa_L./(CLvec*ones(1,length(ts)));
        INa_R_density = INa_R./(CRvec*ones(1,length(ts)));
        
        IK1_ax_density = IK1_ax/Cax;
        IK1_L_density = IK1_L./(CLvec*ones(1,length(ts)));
        IK1_R_density = IK1_R./(CRvec*ones(1,length(ts)));
        
        switch model
            case 'ORd11'
                INaK_ax_density = INaK_ax/Cax;
                INaK_L_density = INaK_L./(CLvec*ones(1,length(ts)));
                INaK_R_density = INaK_R./(CRvec*ones(1,length(ts)));
        end
    end
end
