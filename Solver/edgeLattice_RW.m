function out=edgeLattice_RW(v)

%%% Generate candidate kernel for lattice by random walk.
% Input: edge sets, weight sets for periodic lattices. v must be a squared
% number.
% Output: transition probability matrix
% out(i,j): probability transfer from i to j.
% 

%%%
nsquare = length(v);
out = zeros(nsquare,nsquare);
n = int64(sqrt(nsquare));

rightvector = repmat([ones(1,n-1),0],1,n);
rightvector = rightvector(1:end-1);
out = 0.25*diag(rightvector,1);

makeup_right = repmat([1,zeros(1,n-1)],1,n-1);
makeup_right = [makeup_right,1];
out = out + 0.25*diag(makeup_right,-(n-1));

upvector = ones(1,nsquare-n);
out = out + 0.25*diag(upvector,-n);
makeup_up = ones(1,n);
out = out + 0.25*diag(makeup_up,(n-1)*n);



out = out + out';
% leftvector = repmat([0,ones(1,n-1)],1,n);
% leftvector = leftvector(2:end);
% out = out + 0.25*diag(leftvector,-1);

% out = 0.25*diag(ones(nsquare-1,1),-1);
% out = out + 0.25*diag(ones(nsquare-1,1),1);
% out = out + 0.25*diag(ones(nsquare-n,1),n);
% out = out + 0.25*diag(ones(nsquare-n,1),-n);
% out = out + 0.25*diag([ones(n,1);zeros(nsquare-n-(n-1),1)],n-1);
% out = out + 0.25*diag([ones(n,1);zeros(nsquare-n-(n-1),1)],-(n-1));

% out(1,n) = 0.5;
% out(n,1) = 0.5;

end