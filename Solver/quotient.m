function out=quotient(k)

% generate log-mean of rate of current measure/target measures
% the same function as the score function

N = length(k);
out = zeros(N,N);
for i=1:N
        for j=1:N
            out(i,j) = k(j)/k(i);
        end
end

end