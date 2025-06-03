function out=edgeHyperCube_RW(dimension)

%%% Generate candidate kernel for Hypercube graph by random walk.
% Input: d is the dimension of the HyperCube
% Output: 1/d  *  adjacent matrix of the hyper cube

%%%

N = 2^dimension;
out = zeros(N,N);

for k = 0 : (N-1)

    neighbors = return_neighbor_Hyper_Cube(k, dimension);
    for j = 1 : dimension

        out(k+1, neighbors(j)+1) = 1 / double(dimension);

    end
end

end