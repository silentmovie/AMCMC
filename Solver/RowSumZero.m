function [out]=RowSumZero(A)

%%%
% Input: 
% Output: A Row Sum Zero version of A by replacing its diagonal
% 

%%%
[N,~] = size(A);
for i=1:N
    A(i,i) = -sum(A(i,:))+A(i,i);
end

out = A;

end