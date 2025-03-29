function out=edgeCn_RW(v)

%%%
% Input: edge sets, weight sets. Now this code does not work with input, as
% we deal with Cn graph.
% Output: candidate kernel matrix q_{ij} by random walk.
% out(i,j): probability transfer from i to j.
% 

%%%
n = length(v);
out = zeros(n,n);

out = 0.5*diag(ones(n-1,1),-1);
out = out + 0.5*diag(ones(n-1,1),1);
out(1,n) = 0.5;
out(n,1) = 0.5;

end