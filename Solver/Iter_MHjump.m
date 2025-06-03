function rhohist = Iter_MHjump(pai, rho0, Q, tspan, deltaT, samplesize, seed,mode)

%%% Initialization
%% Data structure:
% 1) pai,rho0 are row vectors
% 2) y(j,:) is the history of k at time j
%    each row of y is a state; each col of y is a time curve of that position

N = length(pai);
TotIt = tspan(2)/deltaT;
rhohist = zeros(TotIt+1, N);
edges = [0.5:1:(N+0.5)];     % bins in the randsample

%% Initial Data
rhohist(1,:) = rho0;         % the first row of rho-history is t=0, the last row is t=tspan(2)

% inpo is array of 1*samplesize, 
% inpo = ransdsample(N,samplesize, true, rho0)';    % entry is the state of one particle 
% current is array of 1*N
% current = histcounts(inpo, N);    % entry is # of particles in each state
current = samplesize / N * ones(1,N);


P = eye(N) + Q*deltaT;       % eq:jump on Page 3, P(i,j) is prob of node i jump to node j.

for j = 2:TotIt+1
    
    if mod(j,1000)==0
        j
    end
    
    % out is matrix of N*N
    out = zeros(N,N);     % out(i,j): # of particles from state i to state j.
    parfor state = 1:N
        %P(state,:) is the transition prob from state to the other probability
        seed_ti = seed + 10000 * j + state;
        stream = RandStream('Threefry', 'Seed', seed_ti);
        tmp = randsample(stream, N, current(state), true, P(state,:))';
        out(state,:) = histcounts(tmp, edges); 
    end
    current = sum(out,1);    % sum along column, current(j) is # of particles in destination state j.
    
    rhohist(j,:) = current/samplesize;
    
end



if strcmp(mode, 'None')
    
   return
   
elseif strcmp(mode, 'Print')
    
   warning('use plot function in TestCases folder')
   return
end

end
