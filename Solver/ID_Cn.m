function [pai,Qrow,minEigQrow] = ID_Cn(seed,N)

%%% Generate pai and Q-matrix for Cn graph by random walk.
rng(seed,"twister")
pai = rand(1,N);                   
pai = pai/sum(pai);
pai = sort(pai,'descend');      % target distribution (not necessary to be 'descend')


Q = AcceptReject_RW(pai);
edge = edgeCn_RW(pai);
Q = Q.*edge;
Qrow = RowSumZero(Q);

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