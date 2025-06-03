function [rhohist, Ham] = Iter_KLjump(pai, rho0, Q, psi0, alphat, tspan, deltaT, samplesize, seed, mode)

N = length(pai);

TotIt = int64(tspan(2)/deltaT);
rhohist = zeros(TotIt+1, N);
Ham = zeros(TotIt+1,1);
edges = [0.5:1:(N+0.5)];

rhohist(1,:) = rho0;
psi = psi0;
k = rho0./pai;

% inpo is array of 1*samplesize, 
% inpo = randsample(N,samplesize, true, rho0)';    % entry is the state of one particle 
% current is array of 1*N
% current = histcounts(inpo, N);    % entry is # of particles in each state
current = samplesize / N * ones(1,N);


Ham(1) = sum(log(rhohist(1,:)).*rhohist(1,:) - log(pai).*rhohist(1,:));
Ham(1) = Ham(1) + sum(0.25*pai*(logmean(k).*psidiffsquare(psi).*Q));



for j = 2:TotIt+1
    j
    out = zeros(N,N);
    W = Q.*(PNpsi(psi))';
    W = diag(pai./rhohist(j-1,:)) * W;
    W = RowSumZero(W);
    W = eye(N) + deltaT * W;
    
    parfor state = 1:N
          % out is matrix of N*N

        seed_ti = seed + 10000 * j + state;
        stream = RandStream('Threefry', 'Seed', seed_ti);
        tmp = randsample(stream, N, current(state), true, P(state,:))';
        out(state,:) = histcounts(tmp, edges);
    end
    
    current = sum(out,1);
    rhohist(j,:) = current/samplesize;
    psi = psi + deltaT*(-alphat*psi -log(k) - 0.5* sum(Q.*psidiffsquare(psi).*partialTheta(k),2)');
    k = rhohist(j,:)./pai;
    if any(k<0) | any(imag(psi)~=0)
        k
        imag(psi)
        break
    end
    Ham(j) = sum(log(rhohist(j,:)).*rhohist(j,:) - log(pai).*rhohist(j,:));
    Ham(j) = Ham(j) + sum(0.25*pai*(logmean(k).*psidiffsquare(psi).*Q));
end

if strcmp(mode, 'None')
    
   return
   
elseif strcmp(mode, 'Print')
    
   warning('use plot function in TestCases folder')
   return
end

end

    