function [Rmat, Cmat, Iind, Nnodes, f_I, iEC, cleft, indices] = ... 
    generate_1D_Mdisc_cleft_EpC(r, L, Ncell, Nint, Mdisc, D, Gb_mat, Gc_array, ...
    IDarea_vec, loc_mat, scaleI, ggap_val)
% generate 1D monodomain matrices, uniform properties
% Rmat, Cmat: Nres x 3, Ncap x 3
% row 1/2: node connections, row 3: value
%
% Iind: Nionic x 2
% row 1/2: node connections
% 
% Nnodes (scalar): number of total nodes

if nargin == 0
    Ncell = 5;
    Nint = 2;   % number of intracellular nodes
    Mdisc = 10; % number of cleft nodes
    
    D = 1; % diffusion coefficient, cm^2/s  - to determing GJ conductance  
    
    Gb_mat = randi(10, 1, Mdisc); % bulk conductances, 1 x Mdisc, mS
    %randomness for testing
    Gb_mat(rand(1,Mdisc)<.5) = 0;
    
    Gc_array = randi(10, Mdisc, Mdisc); % intra-cleft conductances, Mdisc x Mdisc, mS
    %randomness for testing
    Gc_array(rand(Mdisc,Mdisc)<.5) = 0;

    
    IDarea_vec = randi(10, Mdisc, 1); % ID patch surface area, Mdisc x 1, um^2
   
    loc_vec = [.5 0 0];  % localization of the currents at the ID

    tmp = rand(Mdisc,1); tmp = tmp/sum(tmp);
    loc_mat = tmp*loc_vec;
    scaleI = [1 1 1];


end

 % scaling factor for localization of channels at ID
[~,Ncurrents] = size(loc_mat);  % number of currents in ionic model
loc_vec = sum(loc_mat,1);
loc_ax = 1-loc_vec;  % localization on axial membrane


Cm = 1*1e-8;      % membrane capacitance, uF/um^2

Njunctions = Ncell - 1; 
Nnodes = (Nint+2)*Ncell + Njunctions*Mdisc  + 1;  % total nodes, + ground

% cell geometry
% L = 100;        % cell length, um
% r = 11;         % cell radius, um
Aax = 2*pi*r*L; % patch surface area, um^2
Ad = sum(IDarea_vec);    % disc surface area, um^2
Atot = 2*Ad + Aax;  % total surface area, um^2
Apatch = Aax/Nint;

Ctot = Atot*Cm; % total capacitance, uF
Cdisc = IDarea_vec*Cm;  % disc capacitance, Mdisc x 1, uF
Cpatch = Apatch*Cm;   % axial patch capacitance, uF

dx = L/1000;    % mm
Dnorm = D*100/1000;  % mm^2/ms
R = dx^2./Dnorm;  % normalized resistance (ms)
Rgap = R./Ctot;  % k-ohms ( = ms/uF)
ggap = ggap_val.*ones(1,Ncell-1);  % mS, gap junctional conductance, Ncell-1 x 1 vector, between cell i and i+1, for i = 1,...,N-1
% ggap = 1./Rgap.*ones(1,Ncell-1);  % mS, total gap junctional conductance


p_myo = 150*10; % myoplasmic resistivity, k-ohm*um
Rmyo = p_myo*(L/(Nint+1))/(pi*r^2);    % myoplasmic resistance, k-ohm
gmyo = 1/Rmyo;                  % mS, myoplasmic conductance, scalar

iEC = [];
for i = 1:Ncell-1
    for j = 1:Mdisc
        iEC = [iEC (i-1)*(Nint+2+Mdisc)+j+Nint+2];  % indices of extracellular nodes 
    end
end
iEC = [iEC Nnodes];


% resistors
Nres_bulk = nnz(Gb_mat);
Gc_array = triu(Gc_array,1); % only upper diagonals
Nres_cleft = nnz(Gc_array);

Hcleft = zeros(length(iEC)-1, (2*nnz(Gc_array)+nnz(Gb_mat))*(Ncell-1));
Hcleft_all = zeros(4*(length(iEC)-1), 4*(2*nnz(Gc_array)+nnz(Gb_mat))*(Ncell-1));


Nres = Njunctions + Ncell*(Nint+1) + (Nres_cleft + Nres_bulk)*(Ncell-1);

count = 0;
Rmat = zeros(Nres,3);
% myoplasmic resistors
for i = 1:Ncell
    for k = 0:Nint
        count = count + 1;
        Rmat(count, 1) = (i-1)*(Nint+2+Mdisc)+k+1;
        Rmat(count, 2) = (i-1)*(Nint+2+Mdisc)+k+2;
        Rmat(count, 3) = gmyo;
    end
end

% gap junctions
for i = 1:Ncell-1
    count = count + 1;
    Rmat(count, 1) = i*(Nint+2)+(i-1)*Mdisc;
    Rmat(count, 2) = i*(Nint+2+Mdisc)+1; 
    Rmat(count, 3) = ggap(i);
end

% cleft-bulk resistors
for i = 1:Ncell-1
    for j = 1:Mdisc
        if Gb_mat(j)>0
            count = count + 1;
            Rmat(count, 1) = (i-1)*(Nint+2+Mdisc)+j+Nint+2;
            Rmat(count, 2) = Nnodes; 
            Rmat(count, 3) = Gb_mat(j); 
        end
    end
end
% cleft-cleft resistors
cleft_ind_cleft1 = [];
cleft_ind_cleft2 = [];
cleft_g_cleft = [];
for i = 1:Ncell-1
    for j = 1:Mdisc
        for k = j+1:Mdisc
            if Gc_array(j,k)>0
                count = count + 1;
                Rmat(count, 1) = (i-1)*(Nint+2+Mdisc)+j+Nint+2;
                Rmat(count, 2) = (i-1)*(Nint+2+Mdisc)+k+Nint+2; 
                Rmat(count, 3) = Gc_array(j,k); 
                
                cleft_ind_cleft1 = [cleft_ind_cleft1; Rmat(count,1); Rmat(count,2)];
                cleft_ind_cleft2 = [cleft_ind_cleft2; Rmat(count,2); Rmat(count,1)];
                cleft_g_cleft = [cleft_g_cleft; Rmat(count,3); Rmat(count,3)];

                Hcleft(find(iEC == Rmat(count,1)), length(cleft_ind_cleft1)-1) = 1;
                Hcleft(find(iEC == Rmat(count,2)), length(cleft_ind_cleft1)) = 1;

            end
        end
    end
end

% capacitors / ionic currents
Ncap = (Nint + 2*Mdisc)*Ncell - 2*Mdisc + 2; % no cleft on boundaries
Cmat = zeros(Ncap,3);
f_I = zeros(Ncap, Ncurrents);


% axial capacitors
count = 0;  ind_axial = [];
for i = 1:Ncell
    for k = 1:Nint
        count = count + 1;
        Cmat(count, 1) = (i-1)*(Nint+2+Mdisc)+k+1;
        Cmat(count, 2) = Nnodes;
        Cmat(count, 3) = Cpatch;
        ind_axial = [ind_axial count];

        for q = 1:Ncurrents
           f_I(count, q) = scaleI(q)*loc_ax(q)/Nint; 
        end
    end
end

% disc capacitors
ind_disc_pre = [];
ind_disc_post = [];
for i = 1:Ncell
    for k = 1:2 % post-junctional and pre-junctional (left and right sides)
        for j = 1:Mdisc
            count = count + 1;
            
            if k == 1
                ind_disc_post = [ind_disc_post count];
            else
                ind_disc_pre = [ind_disc_pre count];
            end
            if i == 1 && k == 1
                Cmat(count, 1) = 1;
                Cmat(count, 2) = Nnodes;
                Cmat(count, 3) = sum(Cdisc);
                for q = 1:Ncurrents
                    f_I(count, q) = scaleI(q)*loc_vec(q)/2;
                end
                break;
            elseif i == Ncell && k == 2
                Cmat(count, 1) = Nnodes-1;
                Cmat(count, 2) = Nnodes;
                Cmat(count, 3) = sum(Cdisc);
                for q = 1:Ncurrents
                    f_I(count, q) = scaleI(q)*loc_vec(q)/2;
                end
                break;
            else
                if k == 1
                    Cmat(count, 1) = (i-1)*(Nint+2+Mdisc)+1;
                    Cmat(count, 2) = (i-2)*(Nint+2+Mdisc)+j+Nint+2;
                else
                    Cmat(count, 1) = i*(Nint+2)+(i-1)*Mdisc;
                    Cmat(count, 2) = (i-1)*(Nint+2+Mdisc)+j+Nint+2;
                end
            end
            Cmat(count, 3) = Cdisc(j);
            
            for q = 1:Ncurrents % factor of 2 comes from the 2 sides of the cell
                f_I(count, q) = scaleI(q)*loc_mat(j,q)/2;
            end
        end
    end
    
end

% ionic currents - indices for connections
Iind = Cmat(:,1:2);
f_I = f_I'; % Ncurrents x Npatches



% cleft "battery" terms: indices of cleft nodes that couple to the bulk
bulk_ind_cleft1 = Rmat((Rmat(:,2)==Nnodes),1);
bulk_ind_cleft2 = Nnodes*ones(length(bulk_ind_cleft1),1);
bulk_g_cleft = Rmat((Rmat(:,2)==Nnodes),3); % conductances, mS
tmp = [];
for i = 1:length(bulk_ind_cleft1)
    tmp = [tmp find(iEC == bulk_ind_cleft1(i))];
end
Hcleft(tmp,end-length(bulk_ind_cleft1)+1:end) = eye(length(bulk_ind_cleft1));


ind_cleft1 = [cleft_ind_cleft1; bulk_ind_cleft1];
ind_cleft2 = [cleft_ind_cleft2; bulk_ind_cleft2];

Nec = length(iEC)-1; Ni = length(ind_cleft1);
for i = 1:4
    Hcleft_all((i-1)*Nec+1:i*Nec,(i-1)*Ni+1:i*Ni) = Hcleft;
end

cleft.ind_cleft1_phivec = repmat(ind_cleft1,4,1);
cleft.ind_cleft1_Svec = [ind_cleft1; ind_cleft1+Nnodes; ind_cleft1+2*Nnodes; ind_cleft1+3*Nnodes];
cleft.ind_cleft2_phivec = repmat(ind_cleft2,4,1);
cleft.ind_cleft2_Svec = [ind_cleft2; ind_cleft2+Nnodes; ind_cleft2+2*Nnodes; ind_cleft2+3*Nnodes];
cleft.g_cleft = [cleft_g_cleft; bulk_g_cleft];

cleft.Hcleft = sparse(Hcleft_all);


indices.ind_axial = ind_axial;
indices.ind_disc_pre = ind_disc_pre;
indices.ind_disc_post = ind_disc_post;