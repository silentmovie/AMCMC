function [rhohist,psihist, Ham, effSteps,alphathist] = Iter_chi(pai, rho0, Q, psi0, alphat, tspan, deltaT, mode)

%%% Initialization
%% Data structure:
% 1) pai,rho0 are row vectors
% 2) y(j,:) is the history of k at time j
%    each row of y is a state; each col of y is a time curve of that position

N = length(pai);
TotIt = int64(tspan(2)/deltaT);
rhohist = zeros(TotIt+1,N); 
psihist = zeros(TotIt+1,N);
Ham = zeros(TotIt+1,1);

% effSteps(1) and alphathist(1) are dummy.
effSteps = deltaT*ones(TotIt+1,1);
alphathist = alphat*ones(TotIt+1,1);


%% Initial Data           
rhohist(1,:) = rho0;
psihist(1,:) = psi0;
Ham(1) = 0.5*sum((rhohist(1,:)-pai).^2./pai);
Ham(1) = Ham(1) + sum(0.25*pai*(psidiffsquare(psi0).*Q));
% NaN for the dummy value,
effSteps(1) = NaN;
alphathist(1) = NaN;

% rhoCur = rho0;
% psiCur = psi0;
% kCur = rhoCur./psiCur;

%%% Iteration

for  j= 2:(TotIt+1)

    if mod(j,1000)==0
        j
    end
    
    % warm-stary by MH, psihist(j)=-rhohist(j,:)./pai
    if deltaT*double(j) <=0
       
       rhohist(j,:) = rhohist(j-1,:) + deltaT*(rhohist(j-1,:)*Q);
       psihist(j,:) = -rhohist(j,:)./pai;
       Ham(j) = 0.5*sum((rhohist(j,:)-pai).^2./pai);
       Ham(j) = Ham(j) + sum(0.25*pai*(psidiffsquare(psihist(j,:)).*Q));
       alphathist(j) = -1;
       effSteps(j) = deltaT;
       continue
    end


    % Gauss-Seidel iteration
    rhohist(j,:) = rhohist(j-1,:) + deltaT * psihist(j-1,:) * diag(-pai) * Q;
    psihist(j,:) = psihist(j-1,:) + deltaT * (-alphat*psihist(j-1,:) - rhohist(j,:)./pai + 1);
    % psi = psi + deltaT * (-alphat*psi - rhohist(j,:)./pai + 1);
    
    Ham(j) = 0.5*sum((rhohist(j,:)-pai).^2./pai);
    Ham(j) = Ham(j) + sum(0.25*pai*(psidiffsquare(psihist(j,:)).*Q));
end

if strcmp(mode, 'None')
        return
    
elseif strcmp(mode, 'Print')
    
   warning('use plot function in TestCases folder')
   return
end


end