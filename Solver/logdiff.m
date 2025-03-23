function out=logdiff(k)

% generate log-mean of rate of current measure/target measures

N = length(k);
out = zeros(N,N);
for i=1:N
        for j=1:N
            out(i,j) = log(k(i))-log(k(j));
        end
end

end