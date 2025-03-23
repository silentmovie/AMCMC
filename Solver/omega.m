function out=omega(k)

% generate weighted term w.r.t k in the psi equation

N = length(k);
out = zeros(N,N);
for i=1:N
        for j=1:N
            if j == i
                out(i,j) = 0.5;
            elseif abs(k(i)-k(j))/k(j) <= 1e-3
                out(i,j) = 0.5;
            else
            out(i,j) = ( log(k(i)/k(j)) - 1 + k(j)/k(i) ) / (log(k(i)/k(j)))^2;
            end
        end
end

end