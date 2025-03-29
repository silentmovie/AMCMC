function out=AcceptReject_RW(pai, degree)

%%% Accept-Reject matrix with candidate kernel by a random walk
% Input: pai: target distribution in 1*n vector
%        degree: the degress of nodes in 1*n vector
% Output: Accept-Rejection matrix in MH. A_{ij}=min{frac(j)/frac(i), 1}, A_{ii}=0.
%                                        frac(j) = pai(j)/degree(j)

%%%
n = length(pai);


if nargin < 2
   degree = ones(1,n);
end

frac = pai./degree;
out = zeros(n,n);

for i = 1:n
    for j = i+1:n
        out(i,j) = min(frac(j)/frac(i),1);
        out(j,i) = min(frac(i)/frac(j),1);
    end
end

end