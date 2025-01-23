function [P, Q, C] = create_coefficient_matrices(Nnodes, Rmat, Cmat, Iind, dt)
% create coefficient matrices P, Q, C
% inputs: Nnodes, dt, Rmat, Cmat, Iind

% dt = 1e-3; % ms;

[Nres,~] = size(Rmat);
[Ncap,~] = size(Cmat);
[Nionic,~] = size(Iind);

P = sparse(Nnodes, Nnodes);
Q = sparse(Nnodes, Nnodes);
C = sparse(Nnodes, Nionic);



% loop over all resistors
i = Rmat(:,1); j = Rmat(:,2); val = Rmat(:,3);
for k = 1:Nres
    P(i(k),i(k)) = P(i(k),i(k)) + val(k);
    P(i(k),j(k)) = P(i(k),j(k)) - val(k);
    P(j(k),i(k)) = P(j(k),i(k)) - val(k);
    P(j(k),j(k)) = P(j(k),j(k)) + val(k);
end

% loop over all capacitors
i = Cmat(:,1); j = Cmat(:,2); val = Cmat(:,3)/dt;

for k = 1:Ncap
    P(i(k),i(k)) = P(i(k),i(k)) + val(k);
    P(i(k),j(k)) = P(i(k),j(k)) - val(k);
    P(j(k),i(k)) = P(j(k),i(k)) - val(k);
    P(j(k),j(k)) = P(j(k),j(k)) + val(k);
    
    Q(i(k),i(k)) = Q(i(k),i(k)) + val(k);
    Q(i(k),j(k)) = Q(i(k),j(k)) - val(k);
    Q(j(k),i(k)) = Q(j(k),i(k)) - val(k);
    Q(j(k),j(k)) = Q(j(k),j(k)) + val(k);
end

% define last node as ground 
P(Nnodes,:) = 0;
P(Nnodes,Nnodes) = 1;
Q(Nnodes,:) = 0;
Q(Nnodes,Nnodes) = 1;

% loop over all ionic membrane patches
i = Iind(:,1); j = Iind(:,2);
for k = 1:Nionic
    C(i(k),k) = -1;
    C(j(k),k) = 1;
end
C(Nnodes,:) = 0;
