function [pai,Qrow,minEigQrow] = ID_Cn(seed,N)

rng(seed,"twister")
pai = rand(1,N);                   
pai = pai/sum(pai);
pai = sort(pai,'descend');      % target distribution (not necessary to be 'descend')

% some pai used in example.
% 
% pai = [0.456913643618358	0.355387771485700	0.187698584895942];
% rho0 = [0.222092426514290	0.388226387575271	0.389681185910440];

Q = AcceptReject_RW(pai);
edge = edgeCn_RW(pai);
Q = Q.*edge;
Qrow = RowSumZero(Q)

Eigenvalue = eig(Qrow);
minEigQrow = max(Eigenvalue(abs(Eigenvalue)>=1e-3));

% check if detailed balance:
DB = diag(pai)*Qrow;
SymDiff = abs(DB- DB');
if any(SymDiff >=1e-7,"all")
    warning('not detailed balance')
    pause
end

clear DB SymDiff

end