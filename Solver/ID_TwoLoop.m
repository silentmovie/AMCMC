function [pai, Qrow, minEigQrow] = ID_TwoLoop(seed,N)

%%% Generate pai and Q-matrix for twoLoop graph by random walk.
% rng(seed,"twister")
% pai = rand(1,N);                   
% pai = pai/sum(pai);
% degree = 2*ones(1,N);
% degree(floor(N/3))= 3;
% degree(N-floor(N/3))=3;

% manual add target distribution pai and degree matrix
pai = [4/27 4/27 4/27 1/18 1/18 4/27 4/27 4/27];

degree = [2 2 3 2 2 3 2 2];

Q = AcceptReject_RW(pai,degree);
edge = edgeTwoLoop_RW(degree);
Q = Q.*edge;
Qrow = RowSumZero(Q)

Eigenvalue = eig(Qrow);
minEigQrow = max(Eigenvalue(abs(Eigenvalue)>=1e-8));

% check if detailed balance:
DB = diag(pai)*Qrow;
SymDiff = abs(DB- DB');
if any(SymDiff >=1e-7,"all")
    warning('not detailed balance')
    pause
end

clear DB SymDiff

end