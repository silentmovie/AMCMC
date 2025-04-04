function out = psidiffsquare(psi)
% N = length(psi);
% out = zeros(N,N);
% for i=1:N
%     for j=i+1:N
%         out(i,j) = (psi(i)-psi(j)) * (psi(i)-psi(j));
%         out(j,i) = out(i,j);
% 
%     end
% end
% end

psi_diff = psi - psi.';
out = psi_diff.^2;
end