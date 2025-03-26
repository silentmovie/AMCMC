function out=PNpsi(psi)

%%%
% Input: 
% Output: PNpsi(i,j)= Positive of psi(i)-psi(j)

%%%
n = length(psi);

out = zeros(n,n);

for j = 1:n
    for i=1:n
        out(i,j) = max(psi(j)-psi(i),0);
    end
end

end