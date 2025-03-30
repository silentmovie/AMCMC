function [pai,Qrow,minEig] = ID_MGaussian(seed,N)


%%% Generate multi-mode Gaussian on lattices
n = sqrt(N);
x = linspace(0,1,n);
y = linspace(0,1,n);
y = y';
xc1=0.25;
yc1=0.25;
pai = exp(-0.5*((x-xc1).^2+(y-yc1).^2)*20);
xc2=0.75;
yc2=0.75;
pai = pai + exp(-0.5*((x-xc2).^2+(y-yc2).^2)*80);
pai = pai/(sum(sum(pai)));
pai = reshape(pai,1,N);

Q = AcceptReject_RW(pai);
edge = edgeLattice_RW(pai);
Q = Q.*edge;
Qrow = RowSumZero(Q)

Eigenvalue = eig(Qrow);
minEig = max(Eigenvalue(abs(Eigenvalue)>=1e-3));

% check if detailed balance:
DB = diag(pai)*Qrow;
SymDiff = abs(DB- DB');
if any(SymDiff >=1e-7,"all")
    warning('not detailed balance')
    pause
end

clear DB SymDiff

end