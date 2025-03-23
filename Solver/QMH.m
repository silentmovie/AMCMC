function out=QMH(v, degree)

%%%
% Input: any positive vector.
% Output: Q-matrix in MH. Q_{ij}=min{pi_j/pi_i, 1}, Q_{ii}=0.
% 

%%%
n = length(v);

if nargin < 2
   degree = ones(1,n);
end

frac = v./degree;
out = zeros(n,n);

for i = 1:n
    for j = i+1:n
        out(i,j) = min(frac(j)/frac(i),1);
        out(j,i) = min(frac(i)/frac(j),1);
    end
end

end