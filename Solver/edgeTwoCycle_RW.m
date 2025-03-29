function out=edgeTwoCycle_RW(degree)

%%% Generate candidate kernel for TwoCycle Graph by random walk
%
% Input: degree of nodes in 1*n vector,as we deal with TwoCycle graph.
%
% In the first half cycle of size n1, every node is
% connected with two nodes except for the one adjunct node; the adjunct
% node is n-th, from 1 to n-1 in counterclock wise;
% 
% the second half ofcylce of size n2, same;
% the link is single connect from one adjunct node to the other.
% Output: transition probability matrix
% out(i,j): probability transfer from i to j.
% 

%%%
n = length(degree);

weight1 = degree(1);
index = find(degree ~= weight1);
weight2 = degree(index(1));

n1 = index(1);
n2 = n-index(2) + 1;
n3 = n - n1 - n2;

cycle1 = diag(ones(n1-1,1),-1)/weight1;
cycle1 = cycle1 + diag(ones(n1-1,1),1)/weight1;
cycle1(1,n1) = 1/weight1;
cycle1(n1,1) = 1/weight2;
cycle1(n1,n1-1) = 1/weight2;


cycle2 = diag(ones(n2-1,1),-1)/weight1;
cycle2 = cycle2 + diag(ones(n2-1,1),1)/weight1;
cycle2(n2,1) = 1/weight1;
cycle2(1,2) = 1/weight2;
cycle2(1,n2) = 1/weight2;

link = diag(ones(n3-1,1),-1)/weight1;
link = link+ diag(ones(n3-1,1),1)/weight1;
 
out = blkdiag(cycle1, link, cycle2);
out(n1,n1+1)=1/weight2;
out(n1+1,n1)=1/weight1;
out(n1+n3,n1+n3+1)=1/weight1;
out(n1+n3+1,n1+n3)=1/weight2;
out(n1+n3+1,n1+n3+2)=1/weight2;
out(n1+n3+1,n)=1/weight2;
end