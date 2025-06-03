function rhohist=Iter_MH(pai,rho0,Q,tspan,deltaT,mode)

%%% Initialization
%% Data structure:
% 1) pai,rho0 are row vectors
% 2) y(j,:) is the history of k at time j
%    each row of y is a state; each col of y is a time curve of that position

N = length(pai);
TotIt = int64(tspan(2)/deltaT);
rhohist = zeros(TotIt+1,N);            %each row is a state; each col is a time curve of that position


%% Initial Data
rhohist(1,:) = rho0;                   % the first row of rho-history is t=0, the last row is t=tspan(2)                 


for i=2:(TotIt+1)
    % MH-ode defined between detailed-balanced and Q-MH on Page 2
    if mod(i,1000)==0
        i
    end

    rhohist(i,:) = rhohist(i-1,:) + rhohist(i-1,:)*(deltaT*Q);
end

if strcmp(mode, 'None')
    
    return
    
elseif strcmp(mode, 'Print')
    
   warning('use plot function in TestCases folder')
   return
end


end


% 



