function X = constructCell2Patches(Ncell, Mdisc, Y)

Npatches = 2*Mdisc*(Ncell-1) + Ncell; % number of membrane patches
X = nan(Npatches, 1);

ind_Vm = ((1:Ncell)-1)*(2*Mdisc+1)+1;  % indices of axial membrane patch

X(ind_Vm) = Y;



for i = 1:Ncell-1, for j = 1:Mdisc, X((i-1)*(2*Mdisc+1)+j+1) = Y(i); end, end
for i = 2:Ncell, for j = 1:Mdisc, X((i-2)*(2*Mdisc+1)+Mdisc+j+1) = Y(i); end, end

