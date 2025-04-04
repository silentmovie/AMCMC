function out=logmean(k)

% generate log-mean of rate of current measure/target measures

N = length(k);
out = zeros(N,N);
for i=1:N
        for j=i+1:N
            if abs(k(i)-k(j)) < 1e-8
                out(i,j) = k(j);
                out(j,i) = out(i,j);
            else
                out(i,j)= (k(i)-k(j))/log(k(i)/k(j));
                out(j,i)= out(i,j);
            end
        end
end
out = out + diag(k);

end