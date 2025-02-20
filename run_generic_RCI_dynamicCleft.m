clear

% model = 'LR1';
model = 'ORd11';

% tissue = '1D Mdisc cleft EpC';
% tissue = '1D single cleft EpC';
% tissue = '1D Mdisc cleft ID EpC';
tissue = '1D Mdisc cleft ID EpC hetg tissue';

ID_dist = 'chan';
GJ_dist = 'mesh';
scale_gj_loc = 1;
scale_chan_loc = 1;

% time parameters
bcl = 800;  % ms
nbeats = 2;
T = bcl*nbeats;
% T = 20;

% save parameters
save_flag_restart = 0;
save_name_restart = 'baseline_D1_clamp';

save_flag_data = 0;
save_folder = "data/";
save_name_data = save_folder + "test_mesh7.mat";
t_save = [300:10:400];  % ms, time points to save all state variables

% load parameters
load_flag = 0;
load_name = 'Continue_Test';
load_case = 'restart';
load_restart_t = 340; % time (ms) of restart, must be defined if restart using values in "restart" structure

% time step (use different time step between stim and twin)
dt_factor = 1;

dt1 = .01./dt_factor; % ms, dt between stim and twin (0.01 for EpC)
dt2 = .1./dt_factor; % ms, dt between twin and next stim

if scale_chan_loc>=500 || scale_gj_loc>=500
    dt1 = .01./5; % ms, dt between stim and twin (0.01 for EpC)
    dt2 = .1./5; % ms, dt between twin and next stim
end


dtS1 = dt1/5;  % ms, cleft concentration time step 1
dtS2 = dt2/10;   % ms, cleft concentration time step 2

Ns1 = round(dt1/dtS1);  % operator splitting for cleft concentrations
Ns2 = round(dt2/dtS2);  % operator splitting for cleft concentrations
% sampling interval
dt1_samp = dt1*dt_factor*4; % ms 0
dt2_samp = dt2*dt_factor*4; % ms
twin = 50;
trange = [0 T];

% cleft / bulk ionic concentrations;
K_o = 5.4;                  % mM
Na_o = 140;                 % mM
Ca_o = 1.8;                 % mM
clamp_flag = [0; 0; 0; 0]; % Na, K, Ca, A (clamping the cleft), 1 = clamped

ts = get_time_variable(trange, dt1_samp, dt2_samp, twin, bcl);
save_int = 1000;%last x ms to save
ts_save = ts(ts>=(ts(end) - save_int));
Nts = length(ts_save);

% cell geometry
Cm = 1*1e-8;      % membrane capacitance, uF/um^2
L = 100;        % cell length, um
r = 11;         % cell radius, um
Aax = 2*pi*r*L; % patch surface area, um^2
Ad = pi*r^2;    % disc surface area, um^2
Atot = 2*Ad + Aax;  % total surface area, um^2
Ctot = Atot*Cm; % total capacitance, uF

% overall ID localization values [0, 1]
locINa = 0.7; %0.5
locIK1 = 0.2; %0.2
locICa = 0.2; %0.2
locINaK = 0.2; %0.2
locUniform = 2*Ad/Atot;

switch model
    case 'LR1'
        Ncurrents = 6;
        scaleI = ones(1,Ncurrents);
        p.iina = 1; p.iisi = 2; p.iik = 3;
        p.iik1 = 4; p.iikp = 5; p.iib = 6;

        loc_vec = zeros(1, Ncurrents);
        loc_vec(p.iina) = locINa;
        loc_vec(p.iik1) = locIK1;

    case 'ORd11'
        Ncurrents = 14;
        scaleI = ones(1, Ncurrents); % ionic current scaling factors
        p.iina = 1; p.iinal = 2; p.iito = 3;
        p.iical = 4; p.iikr = 5; p.iiks = 6;
        p.iik1 = 7; p.iinaca_i = 8; p.iinaca_ss = 9;
        p.iinak = 10; p.iikb = 11; p.iinab = 12;
        p.iicab = 13; p.iipca = 14;
        
        scaleI(p.iina) = 1;

        loc_vec = zeros(1, Ncurrents); % ID localization vec
        loc_vec(p.iina) = locINa;
        loc_vec(p.iik1) = locIK1;
        loc_vec(p.iinak) = locINaK;
        loc_vec(p.iical) = locICa;

        loc_vec(p.iito) = locUniform;
        loc_vec(p.iikr) = locUniform;
        loc_vec(p.iiks) = locUniform;

    case 'Court98'
        Ncurrents = 13;
        scaleI = ones(1, Ncurrents); % ionic current scaling factors

        p.iina = 1; p.iik1 = 2; p.iito = 3;
        p.iikur = 4; p.iikr = 5; p.iiks = 6;
        p.iibna = 7; p.iibk = 8; p.iibca = 9;
        p.iinak = 10; p.iicap = 11; p.iinaca = 12;
        p.iical = 13;

        loc_vec = zeros(1, Ncurrents); % ID localization vec
        loc_vec(p.iina) = locINa;
        loc_vec(p.iik1) = locIK1;
        loc_vec(p.iinak) = locINaK;
        loc_vec(p.iical) = locICa;

        loc_vec(p.iito) = locUniform;
        loc_vec(p.iikr) = locUniform;
        loc_vec(p.iiks) = locUniform;

end

% cell and tissue parameters
%     FEM_file = 'mesh_data/FEMDATA_V_100Parts_IP16nm_P17nm_Map1_Gjy_Wvy_MJ_0_peri_0_chan_GJ_down_1.5_1.5.mat';
% FEM_file_list =  ['FEMDATA_V_100Parts_IP16nm_P17nm_Map1_Gjy_Wvy_MJ_0_peri_0_chan.mat';...
%                   'FEMDATA_V_100Parts_IP60nm_P60nm_Map1_Gjy_Wvy_MJ_0_peri_0_chan.mat'];
              
FEM_file_list =  {'FEMDATA_12.mat'};
                                                   
mesh_folder = "mesh_data/384/";
load(mesh_folder + FEM_file_list{1}); 
Ncell = 5; % number of cells
Njuncs = Ncell-1;
tissue_legend = zeros(Njuncs,1) + 1; %index that chooses mesh from FEM_file_list; one less node than Ncell
% tissue_legend(8:14) = 2;


switch tissue
    case '1D single cleft EpC'
        Nint = 2;   % number of intracellular nodes
        Mdisc = 1;
        icells = 2;
        D = 1; % cm^2/s
        w = 10e-3;  % cleft width, um
        fVol = 0.1;  % cleft volume scaling factor
        Gc_array = 0;
        [Rmat, Cmat, Iind, Nnodes, f_I, iEC, Vol_cleft, cleft, indices] = ...
            generate_1D_single_cleft_EpC(r, L, Ncell, Nint, D, w, loc_vec, scaleI, fVol);

        Vol_cleft_vec = Vol_cleft*ones(4*(length(iEC)-1),1);

    case '1D Mdisc cleft EpC'
        Nint = 2;   % number of intracellular nodes
        D = 1; % cm^2/s
        fVol = 1;  % cleft volume scaling factor

        icells = 1;

        p_ext = 150*10;  % extracellular resistivity, k-ohm*um
        f_disc = 1; f_bulk = 1; % cleft conductance scaling factors

        Gc_array = f_disc*(FEM_data.cleft_adjacency_matrix)/p_ext;  % mS, Mdisc x Mdisc
        Gb_mat = f_bulk*FEM_data.bulk_adjacency_matrix'/p_ext;  % mS, 1 x Mdisc

        IDarea_vec = FEM_data.partition_surface;  % ID membrane patch surface area, um^2
        
        ggap =  7.35e-04 * D;   % def  7.35e-04  3.6043e-04 1.4168e-04 4.0046e-05 7.9755e-06

        Mdisc = length(Gb_mat);

        %channel localization
        switch ID_dist
            case 'chan'
                loc_mat = zeros(Mdisc, Ncurrents);
                tmp = IDarea_vec; tmp = tmp/sum(tmp);            
                loc_mat(:,:) = loc_vec.*tmp;
                loc_mat(:,p.iina) = loc_vec(p.iina).*FEM_data.Na_area_norm;
                loc_mat(:,p.iinak) = loc_vec(p.iinak).*FEM_data.NKA_area_norm;
                loc_mat(:,p.iik1) = loc_vec(p.iik1).*FEM_data.Kir21_area_norm;
            case 'area'
                tmp = IDarea_vec; tmp = tmp/sum(tmp);
                loc_mat = tmp*loc_vec; % localization proportional to area, Mdisc x Ncurrents matrix
            case 'GJ_single'
                % one GJ plaque closest to center node 
                GJ_area = 100; GJ_adjacent = {[ind_conn find(Gc_array(ind_conn,:))]};
                ggap_array = zeros(Mdisc, 1); % distribute to nodes
                for i = 1:length(GJ_area)
                    ind = GJ_adjacent{i};
                    ggap_array(ind) = ggap_array(ind) + GJ_area(i)/length(ind);
                end
                gj_norm_chan = ggap_array/sum(ggap_array);
                
                loc_mat = zeros(Mdisc, Ncurrents);
                tmp = IDarea_vec; tmp = tmp/sum(tmp);   
                loc_mat(:,:) = loc_vec.*tmp;
                loc_mat(:,p.iina) = loc_vec(p.iina).*gj_norm_chan;
                loc_mat(:,p.iinak) = loc_vec(p.iinak).*gj_norm_chan;
                loc_mat(:,p.iik1) = loc_vec(p.iik1).*gj_norm_chan;
            case 'chan_scaled'
                loc_mat = zeros(Mdisc, Ncurrents);
                tmp = IDarea_vec; tmp = tmp/sum(tmp);            
                loc_mat(:,:) = loc_vec.*tmp;

                Na_chan_norm = FEM_data.Na_area_norm;
                chan_norm_scale = Na_chan_norm + (1-mean(Na_chan_norm));   
                chan_new = (chan_norm_scale.^scale_chan_loc)./(sum(chan_norm_scale.^scale_chan_loc));
                Na_chan_new = chan_new;                       
                loc_mat(:,p.iina) = loc_vec(p.iina).*Na_chan_new;

                NKA_chan_norm = FEM_data.NKA_area_norm;
                chan_norm_scale = NKA_chan_norm + (1-mean(NKA_chan_norm));   
                chan_new = (chan_norm_scale.^scale_chan_loc)./(sum(chan_norm_scale.^scale_chan_loc));
                NKA_chan_new = chan_new;      
                loc_mat(:,p.iinak) = loc_vec(p.iinak).*NKA_chan_new;

                Kir21_chan_norm = FEM_data.Kir21_area_norm;
                chan_norm_scale = Kir21_chan_norm + (1-mean(Kir21_chan_norm));   
                chan_new = (chan_norm_scale.^scale_chan_loc)./(sum(chan_norm_scale.^scale_chan_loc));
                Kir21_chan_new = chan_new;  
                loc_mat(:,p.iik1) = loc_vec(p.iik1).*Kir21_chan_new; 

            case 'chan_scaled_gj_coloc'
                loc_mat = zeros(Mdisc, Ncurrents);
                tmp = IDarea_vec; tmp = tmp/sum(tmp);            
                loc_mat(:,:) = loc_vec.*tmp;

                gj_norm_chan = FEM_data.gj_area_norm;
                gj_norm_scale = gj_norm_chan + (1-mean(gj_norm_chan));   
                gj_new = (gj_norm_scale.^scale_chan_loc)./(sum(gj_norm_scale.^scale_chan_loc));
                gj_norm_chan = gj_new;   

                loc_mat(:,p.iina) = loc_vec(p.iina).*gj_norm_chan;
                loc_mat(:,p.iinak) = loc_vec(p.iinak).*gj_norm_chan;
                loc_mat(:,p.iik1) = loc_vec(p.iik1).*gj_norm_chan;
        end





        [Rmat, Cmat, Iind, Nnodes, f_I, iEC, cleft, indices] = ...
            generate_1D_Mdisc_cleft_EpC(r, L, Ncell, Nint, Mdisc, D, Gb_mat, ...
            Gc_array, IDarea_vec, loc_mat, scaleI,ggap);


        Vol_cleft_vec =  fVol*repmat(FEM_data.partition_volume,4*(Ncell-1),1); % um^3

    case '1D Mdisc cleft ID EpC'
        Nint = 1;   % number of intracellular nodes
        D = 1; % cm^2/s
        flag_compute_ggap = 0;

        baseline_gj_area = 41.66;
        baseline_ggap = 7.35e-04;
        ggap_area_ratio = baseline_ggap./baseline_gj_area;
%             ggap = ggap_area_ratio .* FEM_data.gj_total_area * D
        ggap =  7.35e-04 * D;   % def  7.35e-04  3.6043e-04 1.4168e-04 4.0046e-05 7.9755e-06
        fVol = 1;  % cleft volume scaling factor

        icells = 1;

        p_ext = 150*10;  %*10 extracellular resistivity, k-ohm*um
        f_disc = 1; f_bulk = 1; % cleft conductance scaling factors

        Gc_array = f_disc*(FEM_data.cleft_adjacency_matrix)/p_ext;  % mS, Mdisc x Mdisc
        Gb_mat = f_bulk*FEM_data.bulk_adjacency_matrix'/p_ext;  % mS, 1 x Mdisc

        IDarea_vec = FEM_data.partition_surface;  % ID membrane patch surface area, um^2

        Mdisc = length(Gb_mat);

        rho_ie = 1;  % ratio of intracellular (ID)-to-extracellular (cleft) resistivity

        % GJ area / connection parameters
        [~,ind_conn] = min(sum((FEM_data.partition_centers-mean(FEM_data.partition_centers)).^2,2));

        switch GJ_dist
            case "single"
                % one GJ plaque closest to center node 
                GJ_area = 100; GJ_adjacent = {[ind_conn find(Gc_array(ind_conn,:))]};
                ggap_array = zeros(Mdisc, 1); % distribute to nodes
                for i = 1:length(GJ_area)
                    ind = GJ_adjacent{i};
                    ggap_array(ind) = ggap_array(ind) + GJ_area(i)/length(ind);
                end
                gj_norm = ggap_array/sum(ggap_array);
            case "equal"
               % equal distribution
               GJ_area = ones(1,Mdisc); GJ_adjacent = num2cell(1:Mdisc);
               ggap_array = zeros(Mdisc, 1); % distribute to nodes
                for i = 1:length(GJ_area)
                    ind = GJ_adjacent{i};
                    ggap_array(ind) = ggap_array(ind) + GJ_area(i)/length(ind);
                end
                gj_norm = ggap_array/sum(ggap_array);
            case "mesh"
                gj_norm = FEM_data.gj_area_norm;
            case "mesh_scaled"  %chan - chan_mean +1)^scale - 1 + chan_mean
                gj_norm = FEM_data.gj_area_norm;
                gj_norm_scale = gj_norm + (1-mean(gj_norm));   
                gj_new = (gj_norm_scale.^scale_gj_loc)./(sum(gj_norm_scale.^scale_gj_loc));
                gj_norm = gj_new;               

        end

        %channel localization
        switch ID_dist
            case 'chan'
                loc_mat = zeros(Mdisc, Ncurrents);
                tmp = IDarea_vec; tmp = tmp/sum(tmp);            
                loc_mat(:,:) = loc_vec.*tmp;
                loc_mat(:,p.iina) = loc_vec(p.iina).*FEM_data.Na_area_norm;
                loc_mat(:,p.iinak) = loc_vec(p.iinak).*FEM_data.NKA_area_norm;
                loc_mat(:,p.iik1) = loc_vec(p.iik1).*FEM_data.Kir21_area_norm;
            case 'area'
                tmp = IDarea_vec; tmp = tmp/sum(tmp);
                loc_mat = tmp*loc_vec; % localization proportional to area, Mdisc x Ncurrents matrix
            case 'GJ_single'
                % one GJ plaque closest to center node 
                GJ_area = 100; GJ_adjacent = {[ind_conn find(Gc_array(ind_conn,:))]};
                ggap_array = zeros(Mdisc, 1); % distribute to nodes
                for i = 1:length(GJ_area)
                    ind = GJ_adjacent{i};
                    ggap_array(ind) = ggap_array(ind) + GJ_area(i)/length(ind);
                end
                gj_norm_chan = ggap_array/sum(ggap_array);
                
                loc_mat = zeros(Mdisc, Ncurrents);
                tmp = IDarea_vec; tmp = tmp/sum(tmp);   
                loc_mat(:,:) = loc_vec.*tmp;
                loc_mat(:,p.iina) = loc_vec(p.iina).*gj_norm_chan;
                loc_mat(:,p.iinak) = loc_vec(p.iinak).*gj_norm_chan;
                loc_mat(:,p.iik1) = loc_vec(p.iik1).*gj_norm_chan;
            case 'chan_scaled'
                loc_mat = zeros(Mdisc, Ncurrents);
                tmp = IDarea_vec; tmp = tmp/sum(tmp);            
                loc_mat(:,:) = loc_vec.*tmp;

                Na_chan_norm = FEM_data.Na_area_norm;
                chan_norm_scale = Na_chan_norm + (1-mean(Na_chan_norm));   
                chan_new = (chan_norm_scale.^scale_chan_loc)./(sum(chan_norm_scale.^scale_chan_loc));
                Na_chan_new = chan_new;                       
                loc_mat(:,p.iina) = loc_vec(p.iina).*Na_chan_new;

                NKA_chan_norm = FEM_data.NKA_area_norm;
                chan_norm_scale = NKA_chan_norm + (1-mean(NKA_chan_norm));   
                chan_new = (chan_norm_scale.^scale_chan_loc)./(sum(chan_norm_scale.^scale_chan_loc));
                NKA_chan_new = chan_new;      
                loc_mat(:,p.iinak) = loc_vec(p.iinak).*NKA_chan_new;

                Kir21_chan_norm = FEM_data.Kir21_area_norm;
                chan_norm_scale = Kir21_chan_norm + (1-mean(Kir21_chan_norm));   
                chan_new = (chan_norm_scale.^scale_chan_loc)./(sum(chan_norm_scale.^scale_chan_loc));
                Kir21_chan_new = chan_new;  
                loc_mat(:,p.iik1) = loc_vec(p.iik1).*Kir21_chan_new; 

            case 'chan_scaled_gj_coloc'
                loc_mat = zeros(Mdisc, Ncurrents);
                tmp = IDarea_vec; tmp = tmp/sum(tmp);            
                loc_mat(:,:) = loc_vec.*tmp;

                gj_norm_chan = FEM_data.gj_area_norm;
                gj_norm_scale = gj_norm_chan + (1-mean(gj_norm_chan));   
                gj_new = (gj_norm_scale.^scale_chan_loc)./(sum(gj_norm_scale.^scale_chan_loc));
                gj_norm_chan = gj_new;   

                loc_mat(:,p.iina) = loc_vec(p.iina).*gj_norm_chan;
                loc_mat(:,p.iinak) = loc_vec(p.iinak).*gj_norm_chan;
                loc_mat(:,p.iik1) = loc_vec(p.iik1).*gj_norm_chan;
        end


        [Rmat, Cmat, Iind, Nnodes, f_I, iEC, cleft, indices] = ...
            generate_1D_Mdisc_cleft_ID_EpC(r, L, Ncell, Nint, Mdisc, D, Gb_mat, ...
            Gc_array, IDarea_vec, loc_mat, scaleI, gj_norm, rho_ie,flag_compute_ggap,ggap);

        Vol_cleft_vec =  fVol*repmat(FEM_data.partition_volume,4*(Ncell-1),1); % um^3

    case '1D Mdisc cleft ID EpC hetg tissue'
        
        Nint = 1;   % number of intracellular nodes
        D = 1; % cm^2/s
        flag_compute_ggap = 0;

        baseline_gj_area = 41.66;
        baseline_ggap = 7.35e-04;
        ggap_area_ratio = baseline_ggap./baseline_gj_area;
%             ggap = ggap_area_ratio .* FEM_data.gj_total_area * D
        ggap =  7.35e-04 * D;   % def  7.35e-04  3.6043e-04 1.4168e-04 4.0046e-05 7.9755e-06
        fVol = 1;  % cleft volume scaling factor

        icells = 1;

        p_ext = 150*10;  %*10 extracellular resistivity, k-ohm*um
        f_disc = 1; f_bulk = 1; % cleft conductance scaling factors

        % Gc_array = f_disc*(FEM_data.cleft_adjacency_matrix)/p_ext;  % mS, Mdisc x Mdisc
        % Gb_mat = f_bulk*FEM_data.bulk_adjacency_matrix'/p_ext;  % mS, 1 x Mdisc

        % IDarea_vec = FEM_data.partition_surface;  % ID membrane patch surface area, um^2

        load(mesh_folder + FEM_file_list{tissue_legend(1)}); 

        Mdisc = length(FEM_data.bulk_adjacency_matrix);

        rho_ie = 1;  % ratio of intracellular (ID)-to-extracellular (cleft) resistivity

        % GJ area / connection parameters
        [~,ind_conn] = min(sum((FEM_data.partition_centers-mean(FEM_data.partition_centers)).^2,2));
        
        
        
        
        Gc_array = zeros(Mdisc, Mdisc, Njuncs);
        Gb_mat = zeros(Mdisc, Njuncs);
        IDarea_vec = zeros(Mdisc, 2*Njuncs);
        chan_area_norm_mat = zeros(Mdisc, 2*Njuncs);
        Vol_cleft_vec = [];

        for i = 1:Njuncs
            load(mesh_folder + FEM_file_list{tissue_legend(i)}); 
            Gc_array(:,:,i) = f_disc*(FEM_data.cleft_adjacency_matrix)/p_ext;  % mS, Mdisc x Mdisc
            Gb_mat(:,i) = f_bulk*FEM_data.bulk_adjacency_matrix'/p_ext;  % mS, 1 x Mdisc
            IDarea_vec(:,2*i-1) = FEM_data.partition_surface;  % ID membrane patch surface area, um^2
            IDarea_vec(:,2*i) = FEM_data.partition_surface;  % ID membrane patch surface area, um^2
            Vol_cleft_vec =  [Vol_cleft_vec; fVol*repmat(FEM_data.partition_volume,4,1)]; % um^3
%             gj_norm_list(:,i) = FEM_data.gj_area_norm;
%             chan_area_norm_mat(:,2*i-1) = FEM_data.chan_area_norm;
%             chan_area_norm_mat(:,2*i) = FEM_data.chan_area_norm;
        end

        switch GJ_dist
            case "single"
%                 % one GJ plaque closest to center node 
%                 GJ_area = 100; GJ_adjacent = {[ind_conn find(Gc_array(ind_conn,:))]};
%                 ggap_array = zeros(Mdisc, 1); % distribute to nodes
%                 for i = 1:length(GJ_area)
%                     ind = GJ_adjacent{i};
%                     ggap_array(ind) = ggap_array(ind) + GJ_area(i)/length(ind);
%                 end
%                 gj_norm = ggap_array/sum(ggap_array);
            case "equal"
%                % equal distribution
%                GJ_area = ones(1,Mdisc); GJ_adjacent = num2cell(1:Mdisc);
%                ggap_array = zeros(Mdisc, 1); % distribute to nodes
%                 for i = 1:length(GJ_area)
%                     ind = GJ_adjacent{i};
%                     ggap_array(ind) = ggap_array(ind) + GJ_area(i)/length(ind);
%                 end
%                 gj_norm = ggap_array/sum(ggap_array);

            case "mesh"
                for i = 1:Njuncs
                    load(mesh_folder + FEM_file_list{tissue_legend(i)}); 
                    gj_norm_list(:,i) = FEM_data.gj_area_norm;

                end
                
           case "mesh_scaled"  %chan - chan_mean +1)^scale - 1 + chan_mean
                for i = 1:Njuncs
                    load(mesh_folder + FEM_file_list{tissue_legend(i)}); 
                    gj_norm = FEM_data.gj_area_norm;
                    gj_norm_scale = gj_norm + (1-mean(gj_norm));   
                    gj_new = (gj_norm_scale.^scale_gj_loc)./(sum(gj_norm_scale.^scale_gj_loc));
                    gj_norm = gj_new; 
                    gj_norm_list(:,i) = gj_norm;
                end
        end

        %channel localization
        switch ID_dist
            case 'chan'
                loc_mat = zeros(Mdisc, Ncurrents, 2*Njuncs);
                for i = 1:Njuncs
                    load(mesh_folder + FEM_file_list{tissue_legend(i)});
                    tmp = FEM_data.partition_surface; tmp = tmp/sum(tmp);            
                    
                    %pre junc - def symmetrical
                    loc_mat(:,:,2*i-1) = loc_vec.*tmp;
                    loc_mat(:,p.iina,2*i-1) = loc_vec(p.iina).*FEM_data.Na_area_norm;
                    loc_mat(:,p.iinak,2*i-1) = loc_vec(p.iinak).*FEM_data.NKA_area_norm;
                    loc_mat(:,p.iik1,2*i-1) = loc_vec(p.iik1).*FEM_data.Kir21_area_norm;
                    
                    %post junc
                    loc_mat(:,:,2*i) = loc_vec.*tmp;
                    loc_mat(:,p.iina,2*i) = loc_vec(p.iina).*FEM_data.Na_area_norm;
                    loc_mat(:,p.iinak,2*i) = loc_vec(p.iinak).*FEM_data.NKA_area_norm;
                    loc_mat(:,p.iik1,2*i) = loc_vec(p.iik1).*FEM_data.Kir21_area_norm;
                end
            case 'area'
                tmp = IDarea_vec; tmp = tmp/sum(tmp);
                loc_mat = tmp*loc_vec; % localization proportional to area, Mdisc x Ncurrents matrix
                for i = 1:Njuncs
                    load(mesh_folder + FEM_file_list{tissue_legend(i)});
                    tmp = FEM_data.partition_surface; tmp = tmp/sum(tmp);            
                    %pre junc - def symmetrical
                    loc_mat(:,:,2*i-1) = loc_vec.*tmp;
                    %post junc
                    loc_mat(:,:,2*i) = loc_vec.*tmp;
                end
        end

        [Rmat, Cmat, Iind, Nnodes, f_I, iEC, cleft, indices] = ...
            generate_1D_Mdisc_cleft_ID_EpC_tissue_hetg(r, L, Ncell, Nint, Mdisc, D, Gb_mat, ...
            Gc_array, IDarea_vec, loc_mat, scaleI, gj_norm_list, rho_ie,flag_compute_ggap,ggap);

end

[P1, Q1, C1] = create_coefficient_matrices(Nnodes, Rmat, Cmat, Iind, dt1);
[P2, Q2, C2] = create_coefficient_matrices(Nnodes, Rmat, Cmat, Iind, dt2);
Pinv1 = inv(P1);
Pinv2 = inv(P2);

[Npatches, ~] = size(Iind);

% indices of patches

switch model  
    case 'LR1'
        [p, x0] = InitialConstants_LR91(Atot);
        p.iina = 1; p.iisi = 2; p.iik = 3;
        p.iik1 = 4; p.iikp = 5; p.iib = 6;
        % order is determined by code in fun_name
        % INa, Isi, IK, IK1, IKp, Ib
        Ncurrents = 6;
        scaleI = ones(1,Ncurrents);

        % ionic model-specific parameters
        ionic_fun_name = 'fun_LR1';

        % initial conditions
        Nstate = 8-1;  % number of state variables, excluding Vm, per patch
        p.mLR1 = 1; % flag for modified LR1 model with Ca2+ speedup

        % stimulus parameters
        p.stim_dur = 1;   % ms
        p.stim_amp = .5*80e-8*Atot;    % uA
        p.istim = 1;

        ind_last = (icells-1)*(Nint+2+Mdisc)+Nint+1;
        p.istim = (find(Iind(:,1)<=ind_last & Iind(:,2)==Nnodes));  % indices of axial membrane patch
        p.indstim = ismember((1:Npatches)',p.istim);

        % extracellular ionic concentrations
        p.Na_o = Na_o; p.K_o = K_o; p.Ca_o = Ca_o;

    case 'ORd11'

        % order is determined by code in fun_name
        % INa INaL Ito ICaL IKr IKs IK1 INaCa_i INaCa_ss INaK  IKb INab ICab IpCa
        %          scaleI(2) = 5; scaleI(5) = .15; % generates EADs for bcl = 1000, homogeneous endo, D=1 cable

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

        Nstate = 41-1;  % number of state variables, excluding Vm, per patch

        % stimulus parameters
        p.stim_dur = 2;   % ms
        p.stim_amp = 1*50;    % uA/uF

        ind_last = (icells-1)*(Nint+2+Mdisc)+Nint+1;
        p.istim = (find(Iind(:,1)<=ind_last & Iind(:,2)==Nnodes)); 
        p.indstim = ismember((1:Npatches)',p.istim);

        % cell type
        p.celltype = zeros(Npatches,1); %endo = 0, epi = 1, M = 2
        transmural_flag = 0;

        if transmural_flag
            % transmural wedge: endo -> M -> epi
            % Cells 1:iEndoM -> endo
            % Cells iEndoM+1:iMEpi -> M
            % Cells iMEpi+1:end -> epi
            iEndoM = 60; iMEpi = 105;
            p.celltype(1+iEndoM+(2*iEndoM-1)*Mdisc:iMEpi+(2*iMEpi-1)*Mdisc) = 2;
            p.celltype(1+iMEpi+(2*iMEpi-1)*Mdisc:end) = 1;
        end

    case 'Court98'
        % order is determined by code in fun_name
        % INa IK1 Ito IKur IKr IKs IBNa IBK IBCa INaK ICaP INaCa ICaL



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

        Nstate = 21-1;  % number of state variables, including Vm, per patch
        % stimulus parameters
        p.stim_dur = 1;   % ms
        p.stim_amp = -60;    % uA/uF

        ind_last = (icells-1)*(Nint+2+Mdisc)+Nint+1;
        p.istim = (find(Iind(:,1)<=ind_last & Iind(:,2)==Nnodes));   % indices of axial membrane patches
        p.indstim = ismember((1:Npatches)',p.istim);       
end

% cleft / bulk ion concentrations
A_o = Na_o + K_o + 2*Ca_o;  % anion A- concentration, mM

p.bcl = bcl;
p.Npatches = Npatches;
p.L = L; p.r = r;
p.Ctot = Atot*Cm;   % total cell capacitance, uF

Ncleft_comp = length(iEC)-1;

% initalize / load initial conditions
phi_mat = nan(Nnodes, Nts);
G_mat = nan(Nstate*Npatches, Nts);
S_mat = nan(4*Nnodes, Nts); % Na, K, Ca, A
I_all = nan(Npatches*Ncurrents, Nts);

if save_flag_restart
   phi_save = nan(Nnodes, length(t_save));
   G_save = nan(Nstate*Npatches, length(t_save));
   S_save = nan(4*Nnodes, length(t_save));
   I_save = nan(Npatches*Ncurrents, length(t_save));
   count_save = 1;
end

if ~load_flag
    phi0 = x0(1)*ones(Nnodes,1); % intracellular nodes
    phi0(iEC) = 0; % extracellular nodes

    Scleft = nan(4*Nnodes, 1); zvec = nan(4*Nnodes, 1);
    Scleft(1:Nnodes) = Na_o; zvec(1:Nnodes) = 1;
    Scleft(Nnodes+1:2*Nnodes) = K_o; zvec(Nnodes+1:2*Nnodes) = 1;
    Scleft(2*Nnodes+1:3*Nnodes) = Ca_o; zvec(2*Nnodes+1:3*Nnodes) = 2;
    Scleft(3*Nnodes+1:4*Nnodes) = A_o; zvec(3*Nnodes+1:4*Nnodes) = -1;
    S_mat(:,1) = Scleft;

    g0 = nan(Nstate*Npatches,1);
    for i = 1:Nstate
        g0(Npatches*(i-1)+1:i*Npatches) = x0(i+1);
    end
else
    load(load_name);
    switch load_case
        case 'final'
            phi0 = final.phi;
            g0 = final.G;
            Scleft = final.S;
            ts = ts + final.t;
        case 'restart'
           ts = ts + load_restart_t;
           [~,i] = min(abs(restart.t - load_restartq_t));
           phi0 = restart.phi(:,i);
           g0 = restart.G(:,i);
           Scleft = restart.S(:,i);
    end

end

phi_mat(:,1) = phi0;
G_mat(:,1) = g0;

% constants
F = 96.5;                   % Faraday constant, coulombs/mmol
R = 8.314;                  % gas constant, J/K
Temp = 273+37;                 % absolute temperature, K
RTF=(R*Temp/F);                % mV

% cleft concentration clamping
clamp_vec = ones(4*Nnodes, 1);
for i = 1:4
    clamp_vec(iEC(1:end-1) + (i-1)*Nnodes) = clamp_flag(i);
end

% "disc" currents: indices of membrane patches that couple to corresponding
% cleft nodes/compartments (by design, should always be 2 per cleft node,
% pre- and post-junctional membrane pathces)
ind_disc = nan(2, length(iEC)-1);
for i = 1:length(iEC)-1
    ind_disc(:,i) = find(Iind(:,2)==iEC(i));
end
q_disc = [iEC(1:end-1) iEC(1:end-1)+Nnodes iEC(1:end-1)+2*Nnodes iEC(1:end-1)+3*Nnodes];
ind_cleft1_phivec = cleft.ind_cleft1_phivec; ind_cleft2_phivec = cleft.ind_cleft2_phivec;
ind_cleft1_Svec = cleft.ind_cleft1_Svec; ind_cleft2_Svec = cleft.ind_cleft2_Svec;
g_cleft_vec = repmat(cleft.g_cleft, 4, 1)/sum(~clamp_flag);
if ~sum(~clamp_flag)
    g_cleft_vec = 0;
end
Ng = length(cleft.g_cleft); zall = ones(4*Ng, 1); zall(2*Ng+1:3*Ng) = 2; zall(3*Ng+1:4*Ng) = -1;
Hcleft = cleft.Hcleft;

H = zeros(Nnodes, length(iEC)-1);
for i = 1:length(iEC)-1
    H(iEC(i),i) = 1;
end
H = sparse(H);

ionic_fun = str2func(['@(t,x,p,S) ',ionic_fun_name,'(t,x,p,S)']);
p.f_I = f_I;

count = 1; phi_i = phi_mat(:,1); G_i = G_mat(:,1);
% collect Vm
Vm = phi_i(Iind(:,1)) - phi_i(Iind(:,2));
Sp = [Scleft(Iind(:,2)); Scleft(Iind(:,2)+Nnodes); Scleft(Iind(:,2)+2*Nnodes)];
[~, ~, ~, ~, I_new] = ionic_fun(0, [Vm; G_i], p, Sp);
I_all(:,1) = I_new;

beat_num = ones(Npatches,1);
tup = nan(Npatches,1);
trepol = nan(Npatches,1);
Vm_old = Vm; Vthresh = -60; % mV

ti = 0;  % initialize time
tic
while ti < T

    % display
    if ~mod(ti,50)
        disp(ti);
    end

    if mod(ti,bcl)<twin
        dt = dt1; dt_samp = dt1_samp; Ns = Ns1; Q = Q1; C = C1; Pinv = Pinv1;
    else
        dt = dt2; dt_samp = dt2_samp; Ns = Ns2; Q = Q2; C = C2; Pinv = Pinv2;
    end
    p.dt = dt;

    % collect Vm
    Vm = phi_i(Iind(:,1)) - phi_i(Iind(:,2));

    % collect S (cleft ionic concentration)
    Sp = [Scleft(Iind(:,2)); Scleft(Iind(:,2)+Nnodes); Scleft(Iind(:,2)+2*Nnodes)];
    % calculate ionic currents / update gating variables
    [G_new, Iion, Ivec, ~, I_new] = ionic_fun(ti, [Vm; G_i], p, Sp);


    % update cleft currents
    % disc currents
    Idisc = [sum(Ivec([ind_disc, ind_disc + Npatches, ind_disc + 2*Npatches])), zeros(1, length(iEC)-1)]';

    %  cleft-cleft and cleft-bulk currents
    for i = 1:Ns
        Erev = RTF./zall.*(log(Scleft(ind_cleft2_Svec)./Scleft(ind_cleft1_Svec)));
        Ibulk_term = Hcleft*(g_cleft_vec.*Erev);
        Ibulk = Hcleft*(g_cleft_vec.*(phi_i(ind_cleft1_phivec) - phi_i(ind_cleft2_phivec) - Erev));
        dS = ~clamp_vec(q_disc)*dt/Ns.*((Idisc-Ibulk)*1e6./(zvec(q_disc)*F.*Vol_cleft_vec));

        Scleft(q_disc) = Scleft(q_disc) + dS;
    end

    Hbulk = sum(reshape(Ibulk_term, length(iEC)-1,4),2);
    % update voltages
    phi_new = Pinv*(Q*phi_i + C*Iion - H*Hbulk);

      % calculate activation/repolarization times
    [tup, trepol, beat_num] = update_tup_repol(ti, dt, Vm, Vm_old, Vthresh, ...
        tup, trepol, beat_num);

    Vm_old = Vm;
    ti = round(ti + dt, 5);
    
    if any(isnan(Vm))
        disp("complex Vm")
        break
    end

    if ~mod(ti, dt_samp) && ti>(ts(end) - save_int)
        phi_mat(:,count+1) = phi_new;
        G_mat(:,count+1) = G_new;
        S_mat(:,count+1) = Scleft;
        I_all(:,count) = I_new;
        count = count + 1;
    end
    phi_i = phi_new;
    G_i = G_new;

    if save_flag_restart
       if ismember(ti, t_save)
          phi_save(:,count_save) = phi_new;
          G_save(:,count_save) = G_new;
          S_save(:,count_save) = Scleft;
          I_save(:,count_save) = I_new;
          count_save = count_save + 1;
       end
    end
end
toc
%% plotting 
switch tissue

    case '1D single cleft EpC'
        icleft = iEC(1:end-1);
        iintra = setdiff(1:Nnodes-1,icleft);
        phi_i = phi_mat(iintra,:);
        phi_cleft = phi_mat(icleft,:);

        % collect / sort Vm, ionic currents
        Vm = phi_mat(Iind(:,1),:) - phi_mat(Iind(:,2),:);

        [~,ind] = sort(Iind(:,1));
        Iind_cable = Iind(ind,:);
        Vm_cable = Vm(ind,:);
        tup = tup(ind,:); trepol = trepol(ind,:);


        INa_all = I_all(p.iina:Ncurrents:end,:);
        INa_all = INa_all(ind,:); % sort to match Vm order

        [ind_axial, ~] = ind2sub([length(ind) length(indices.ind_axial)], find(ind == indices.ind_axial));
        [ind_disc_pre, ~] = ind2sub([length(ind) length(indices.ind_disc_pre)], find(ind == indices.ind_disc_pre));
        [ind_disc_post, ~] = ind2sub([length(ind) length(indices.ind_disc_post)], find(ind == indices.ind_disc_post));
        INa_axial = INa_all(ind_axial,:);
        INa_disc_pre = INa_all(ind_disc_pre,:);
        INa_disc_post = INa_all(ind_disc_post,:);

        subplot(3,1,1); h = pcolor(ts, 1:Nnodes-1, phi_mat(1:end-1,:));
        set(gca,'ydir','reverse');
        h.LineStyle = 'none';
        subplot(3,1,2); plot(ts, phi_i); hold on;
        subplot(3,1,3); plot(ts, phi_cleft); hold on;

        tup_axial = tup(ind_axial(1:Nint:end),:);
        i1 = round(.25*Ncell); i2 = round(.75*Ncell);
        cv_est = 100*(i2-i1)*(p.L/1000)./(tup_axial(i2,:)-tup_axial(i1,:));  % mm/ms = m/s, 100*m/s = cm/s

    case {'1D Mdisc cleft EpC','1D Mdisc cleft ID EpC', '1D Mdisc cleft ID EpC hetg tissue'}
        %%
        s = 'k-';
        icleft = iEC(1:end-1);
        iintra = setdiff(1:Nnodes-1,icleft);
        phi_i = phi_mat(iintra,:);
        phi_cleft = phi_mat(icleft,:);

        % collect Vm / ionic currents
        Vm = phi_mat(Iind(:,1),:) - phi_mat(Iind(:,2),:);

        [~,ind] = sort(Iind(:,1));
        Iind_cable = Iind(ind,:);
        Vm_cable = Vm(ind,:);
        tup = tup(ind,:); trepol = trepol(ind,:);

        INa_all = I_all(p.iina:Ncurrents:end,:);
        INa_all = INa_all(ind,:); % sort to match Vm order

        % sort axial / pre-/post-junctional membrane indices to match Vm order
        [ind_axial, ~] = ind2sub([length(ind) length(indices.ind_axial)], find(ind == indices.ind_axial));
        [ind_disc_pre, ~] = ind2sub([length(ind) length(indices.ind_disc_pre)], find(ind == indices.ind_disc_pre));
        [ind_disc_post, ~] = ind2sub([length(ind) length(indices.ind_disc_post)], find(ind == indices.ind_disc_post));
        INa_axial = INa_all(ind_axial,:);
        INa_disc_pre = INa_all(ind_disc_pre,:);
        INa_disc_post = INa_all(ind_disc_post,:);
        % Note there is only 1 post-junctional membrane for cell 1 and 
        % 1 pre-junctional membrane for cell N for all values of Mdisc, 
        % which is coupled directly to the bulk
        % To exclude this junction, remove first / last membrane for post / pre

        % indices for a specific junction (pre and post)
        Njunc = 3;
        ind_post = ind_disc_post((Njunc-1)*Mdisc+2:Njunc*Mdisc+1);
        Vm_post = Vm_cable(ind_post,:); INa_post = INa_all(ind_post,:);
        ind_pre = ind_disc_pre((Njunc-1)*Mdisc+1:Njunc*Mdisc);
        Vm_pre = Vm_cable(ind_pre,:); INa_pre = INa_all(ind_pre,:);

        Na_cleft_all = S_mat(iEC,:);
        Na_cleft = Na_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:);       
        K_cleft_all = S_mat(iEC+Nnodes,:);
        K_cleft = K_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:);     
        Ca_cleft_all = S_mat(iEC+2*Nnodes,:);
        Ca_cleft = Ca_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:); 
        A_cleft_all = S_mat(iEC+3*Nnodes,:);
        A_cleft = A_cleft_all((Njunc-1)*Mdisc+1:Njunc*Mdisc,:);
        
        
        

%             figure(100); subplot(2,1,1);F
%             plot(ts, Vm_post,'k',ts, Vm_post(ind_conn,:),'r--'); title('single GJ');
% 



        tup_axial = tup(ind_axial(1:Nint:end),:);
        i1 = round(.25*Ncell); i2 = round(.75*Ncell);
        cv_est = 100*(i2-i1)*(p.L/1000)./(tup_axial(i2,:)-tup_axial(i1,:));  % mm/ms = m/s, 100*m/s = cm/s
        
        
        
        phi_i = phi_mat(iintra,:);
        phi_axial = phi_i(ind_axial,:); 
% 
%         
        figure
        plot(ts_save,Ca_cleft(:,1:end),'color',[1 0 0 0.5])
        ylim([1,3])
%         
        figure
        plot(ts_save,Na_cleft(:,1:end),'color',[0 0 1 0.5])
        ylim([120,150])
        

        for i = 1:25
            phi = phi_axial(i,:);
            phi(phi>-70) = 1;
            phi(phi<=-70) = 0;
            [~,peak_time] = findpeaks(phi);
            
            peak_list(i) = ts_save(peak_time(end));
        end
        
        figure
        hold on
        for i = 1:25
            plot(ts_save,phi_axial(i,1:end))
        end
        
        
        CV = 0.01./(diff(peak_list)./1000);
%         
%         figure
%         plot(ts_save,K_cleft(:,1:end-1),'color',[0 1 0 0.5])
%         ylim([4.5,8])
%         
%         figure
%         plot(ts_save,A_cleft(:,1:end-1),'color',[0 0 0 0.5])
%         ylim([140,150])
% 
        figure
        plot(ts_save,phi_axial(5,1:end))
% 
%         figure;
%         h = pcolor(ts, 1:length(iintra), phi_mat(iintra,:));
%         set(gca,'ydir','reverse');
%         h.LineStyle = 'none';
%         set(gca,'ytick',1:length(iintra), 'yticklabel',iintra);
%         colorbar
%         caxis([-80,0])
% 
        % subplot(2,2,3); plot(ts_save, phi_i,s); hold on;
% 
%         subplot(2,2,2);
%         h = pcolor(ts, 1:length(icleft), phi_mat(icleft,:));
%         set(gca,'ydir','reverse');
%         h.LineStyle = 'none';
%         set(gca,'ytick',1:20:length(icleft), 'yticklabel',icleft(1:20:end));
% 
%         subplot(2,2,4); plot(ts, phi_cleft, s); hold on;
        
%         figure
%         plot(INa_all(ind_pre,:)'.*1e6)
%         ylim([-140,20])
%%

end

%% analysis
  
% list_higher_post = [];
% list_higher_pre = [];
% for i = 1:100
%     
%     if min(INa_pre(i,1:end)) > min(INa_post(i,1:end))
%         list_higher_post = [list_higher_post, i];
%     else
%         list_higher_pre = [list_higher_pre, i];
%     end
% end
% 
% figure
% plot3(FEM_data.partition_centers(list_higher_post, 1), ...
%       FEM_data.partition_centers(list_higher_post, 2), ...
%       FEM_data.partition_centers(list_higher_post, 3), '.','color','red', 'markersize',10)
% hold on
% plot3(FEM_data.partition_centers(list_higher_pre, 1), ...
%       FEM_data.partition_centers(list_higher_pre, 2), ...
%       FEM_data.partition_centers(list_higher_pre, 3), '.','color','blue', 'markersize',10)
% 
%   
%   
% figure
% plot(ts_save,INa_pre(list_higher_post,1:end),'color',[0 0 1 0.5])
% hold on
% plot(ts_save,INa_post(list_higher_post,1:end),'color',[0 1 1 0.5])
% xlim([400.8,401.6])
% 
% figure
% plot(ts_save,Na_cleft(list_higher_post,1:end-1),'color',[0 0 1 0.5])
% hold on
% plot(ts_save,Na_cleft(list_higher_pre,1:end-1),'color',[0 1 1 0.5])
% 
% xlim([400.8,405.6])


    
%% saving
if save_flag_restart
   % add all variables and/or parameters to save,
   % For restarting (continuing):
   % "final" will include final states
   % "restart" will include all states for times in "t_save"
   
   final.phi = phi_new;
   final.G = G_new;
   final.S = Scleft;
   final.I = I_new;
   final.t = ts(end);

   restart.t = t_save;
   restart.phi = phi_save;
   restart.G = G_save;
   restart.S = S_save;
   restart.I = I_save;

   save(save_name_restart,'restart', 'final');

end

if save_flag_data   
   p.loc_vec = loc_vec;
   %tup and trepol are already indexed
   save(save_name_data,'p','iEC','Nnodes','Ncell','Ncurrents','indices','Mdisc',...
       'phi_mat','Iind','S_mat','I_all','ts','model','FEM_file_list',...
       'tissue_legend','tup','trepol','cv_est','ts_save','Nint');
end







