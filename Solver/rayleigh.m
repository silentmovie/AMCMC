function out=rayleigh(Omega,pai)

N = length(pai);

%Create a (N*(N-1)) ONB that spans the space orthogonal to (1,1,...1)

V = [ones(N-1,1),-eye(N-1)];
invPai = diag(1./pai);
V = V';
V = V/sqrt(2);
R= chol(V'*Omega*V);
invR = inv(R);
Mat4Eig = invR'*(V'*Omega*invPai*Omega*invPai*Omega*V)*invR;
Eigenvalue = eig(Mat4Eig);
Eigenvalue = sort(Eigenvalue,'ascend');
out = Eigenvalue;
end


% N = length(pai);
% 
% %Create a (N*(N-1)) ONB that spans the space orthogonal to (1,1,...1)
% 
% V = [ones(N-1,1),-eye(N-1)];
% invPai = diag(1./pai);
% V = V';
% V = V/sqrt(2);
% R= chol(V'*(-Omega)*V);
% invR = inv(R);
% Mat4Eig = invR'*(V'*Omega*invPai*(-Omega)*invPai*Omega*V)*invR;
% Eigenvalue = eig(Mat4Eig);
% Eigenvalue = sort(Eigenvalue,'ascend');
% out = Eigenvalue;
% end


