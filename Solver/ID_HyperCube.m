function [pai,Qrow,minEigQrow] = ID_HyperCube(seed,N)


rng(seed,"twister")

d = int64(log2(N))
% pai is unnormalized!
pai = ones(1,N);                   
pai(1) = 16;
pai(N) = 16;

Q = AcceptReject_RW(pai);
edge = edgeHyperCube_RW(d);
Q = Q.*edge;
Qrow = RowSumZero(Q);

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