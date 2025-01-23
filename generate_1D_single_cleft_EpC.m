function [Rmat, Cmat, Iind, Nnodes, f_I, iEC, Vol_cleft, cleft, indices] = ... 
    generate_1D_single_cleft_EpC(r, L, Ncell, Nint, D, w, loc_vec, scaleI, fVol)

% generate 1D monodomain matrices, uniform properties
% Rmat, Cmat: Nres x 3, Ncap x 3
% row 1/2: node connections, row 3: value
%
% Iind: Nionic x 2
% row 1/2: node connections
% 
% Nnodes (scalar): number of total nodes
% 
% iEC: vector of extracellular nodes
%


if nargin == 0
    Ncell = 5;
    Nint = 2;   % number of intracellular nodes
    
    D = 1; % cm^2/s    
    w = 10e-3;  % cleft width, um
    loc_vec = [.9 0 0];  % localization of the currents at the ID
    scaleI = [1 1 1];
end
if nargin < 7
    fVol = 1;
end

Cm = 1*1e-8;      % membrane capacitance, uF/um^2
Ncurrents = length(loc_vec);

Njunctions = Ncell - 1; 
Nnodes = (Nint+2)*Ncell + Njunctions  + 1;  % total nodes, + ground

% cell geometry
% L = 100;        % cell length, um
% r = 11;         % cell radius, um
Aax = 2*pi*r*L; % patch surface area, um^2
Ad = pi*r^2;    % disc surface area, um^2
Atot = 2*Ad + Aax;  % total surface area, um^2
Apatch = Aax/Nint;
Vol_cleft = fVol*Ad*w; % cleft volume, um^3

Ctot = Atot*Cm; % total capacitance, uF
Cdisc = Ad*Cm;  % disc capacitance, uF
Cpatch = Apatch*Cm;   % axial patch capacitance, uF

dx = L/1000;    % mm
Dnorm = D*100/1000;  % mm^2/ms
R = dx^2./Dnorm;  % normalized resistance (ms)
Rgap = R./Ctot;  % k-ohms ( = ms/uF)
ggap = 1./Rgap*ones(1,Ncell-1);  % mS, gap junctional conductance, Ncell-1 x 1 vector, between cell i and i+1, for i = 1,...,N-1

p_myo = 150*10; % myoplasmic resistivity, k-ohm*um
Rmyo = p_myo*(L/(Nint+1))/(pi*r^2);    % myoplasmic resistance, k-ohm
gmyo = 1/Rmyo;                  % mS, myoplasmic conductance, scalar

p_ext = 150*10;  % extracellular resistivity, k-ohm*um
Rradial = p_ext/(8*pi*w);   % radial cleft resistance, k-ohm
gbulk = 1/Rradial; % mS, radial cleft conductance, 


iEC = (Nint+3):(Nint+3):Nnodes;  % indices of extracellular nodes 

% resistors
Nres = 2*Njunctions + Ncell*(Nint+1);

count = 0;
Rmat = zeros(Nres,3);
% myoplasmic resistors
for i = 1:Ncell
    for j = 1:Nint+1
        count = count + 1;
        Rmat(count, 1) = j + (i-1)*(Nint+3);
        Rmat(count, 2) = j + 1 + (i-1)*(Nint+3);
        Rmat(count, 3) = gmyo;
    end
end
% cleft resistors
for i = 1:Ncell-1
    count = count + 1;
    Rmat(count, 1) = (Nint+3)*i;
    Rmat(count, 2) = Nnodes; 
    Rmat(count, 3) = gbulk; 
end
% gap junctions
for i = 1:Ncell-1
    count = count + 1;
    Rmat(count, 1) = (Nint+3)*i-1;
    Rmat(count, 2) = (Nint+3)*i+1; 
    Rmat(count, 3) = ggap(i);
end

% capacitors / ionic currents
Ncap = (Nint+2)*Ncell;
Cmat = zeros(Ncap,3);
f_I = zeros(Ncap, Ncurrents);


% axial capacitors
count = 0; ind_axial = [];
for i = 1:Ncell
    for j = 1:Nint
        count = count + 1;
        Cmat(count, 1) = j + (i-1)*(Nint+3) + 1;
        Cmat(count, 2) = Nnodes;
        Cmat(count, 3) = Cpatch;
        ind_axial = [ind_axial count];
        
        for q = 1:Ncurrents
           f_I(count, q) = scaleI(q)*(1-loc_vec(q))/Nint; 
        end
    end
end

% disc capacitors
ind_disc_pre = [];
ind_disc_post = [];
for i = 1:Ncell
    for j = 1:2
        count = count + 1;
        
        if j == 1
            ind_disc_post = [ind_disc_post count];
        else
            ind_disc_pre = [ind_disc_pre count];
        end
        if i == 1 && j == 1
            Cmat(count, 1) = 1;
            Cmat(count, 2) = Nnodes;
        elseif i == Ncell && j == 2
           Cmat(count, 1) = Nnodes-1;
           Cmat(count, 2) = Nnodes;
       else
           if j == 1
               Cmat(count, 1) = (i-1)*(Nint+3) + 1;
               Cmat(count, 2) = (i-1)*(Nint+3);
           else
               Cmat(count, 1) = i*(Nint+3) - 1;
               Cmat(count, 2) = i*(Nint+3);
           end
       end
       Cmat(count, 3) = Cdisc;

       
        for q = 1:Ncurrents
           f_I(count, q) = scaleI(q)*loc_vec(q)/2; 
        end
    end
   
end

% ionic currents - indices for connections
Iind = Cmat(:,1:2);
f_I = f_I'; % Ncurrents x Npatches

indices.ind_axial = ind_axial;
indices.ind_disc_pre = ind_disc_pre;
indices.ind_disc_post = ind_disc_post;

% cleft "battery" terms: indices of cleft nodes that couple to the bulk
ind_cleft1 = Rmat((Rmat(:,2)==Nnodes),1);
cleft.ind_cleft1_phivec = repmat(ind_cleft1,4,1);
cleft.ind_cleft1_Svec = [ind_cleft1; ind_cleft1+Nnodes; ind_cleft1+2*Nnodes; ind_cleft1+3*Nnodes];
cleft.ind_cleft2_phivec = repmat(Nnodes*ones(length(iEC)-1,1),4,1);
cleft.ind_cleft2_Svec = [Nnodes*ones(length(ind_cleft1),1); 2*Nnodes*ones(length(ind_cleft1),1); 3*Nnodes*ones(length(ind_cleft1),1); 4*Nnodes*ones(length(ind_cleft1),1)];
cleft.g_cleft = Rmat((Rmat(:,2)==Nnodes),3); % conductances, mS
cleft.Hcleft = sparse(eye(4*length(ind_cleft1)));
