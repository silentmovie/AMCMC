function rhohist = Iter_MHjump(pai, rho0, Q, tspan, deltaT, samplesize, mode)

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
inpo = randsample(N,samplesize, true, rho0)';    % entry is the state of one particle 
% current is array of 1*N
current = histcounts(inpo, N);    % entry is # of particles in each state



P = eye(N) + Q*deltaT;       % eq:jump on Page 3, P(i,j) is prob of node i jump to node j.

for j = 2:TotIt+1
    
    if mod(j,1000)==0
        j
    end
    
    % out is matrix of N*N
    out = zeros(N,N);     % out(i,j): # of particles from state i to state j.
    parfor state = 1:N
        %P(state,:) is the transition prob from state to the other probability
       
        tmp = randsample(N, current(state), true, P(state,:))';
        out(state,:) = histcounts(tmp, edges); 
    end
    current = sum(out,1);    % sum along column, current(j) is # of particles in destination state j.
    
    rhohist(j,:) = current/samplesize;
    
end


% error = sqrt(sum((rhohist(:,1:N)-pai).^2,2));     % cal row norm of the difference
% logError = log10(error);

% %%%%%%%%%%%
% size(logError)

if strcmp(mode, 'None')
    
   return
   
elseif strcmp(mode, 'Print')
    
    cc = hsv(3*N);
    t = [tspan(1):deltaT:tspan(2)]';
    
    pltstep = 1;         % plt every pltstep iterations
    figure('Renderer', 'painters', 'Position', [10 10 1200 1200])
    
    subplot(2,1,1)
%     hold on;
%     for state = 1:N
%         plot((0:TotIt),pai(state), 'color', cc(3*(state-1)+1,:),'marker','.')
%         plot((0:pltstep:TotIt), rhohist(1:pltstep:end, state), 'color', cc(3*(state-1)+1,:),'marker','*');
%     end
    
%     xlim([0,tspan(2)])
    ylim([-0.1,1])
    % title('3-point convergence via M-H')
%     title({[num2str(N),'-point convergence time stepsize ', num2str(deltaT, '%0.2e')]
%     ['largest negative eigen=', num2str(minEig, '%0.2e')]
%     });
    hold off;

    error = sqrt(sum((rhohist(:,1:N)-pai).^2,2));     % cal row norm of the difference
    logError = log10(error);

    %%%%%%%%%%%
    size(logError)

    subplot(2,1,2)
    % plot(log(t(1:(end-1)/10000:end)),logError(1:(end-1)/10000:end))
    hold on;
    plot(t,logError)
    plot(t,-3 + 0.0 * t)
    xlim(tspan)
    ylim([-4,0])
    % title('M-H: log-log plot of error');
    title({['MH: t-log(error) plot']
        ['timestep=', num2str(deltaT, '%0.2e')]
        });
    halftime = tspan(2)/2;
    P = polyfit(t(t>halftime),logError(t>halftime),1);
    plot(t(t>halftime), P(1)*t(t>halftime) + P(2), 'r', 'LineWidth', 2)
    % Define the formula using LaTeX syntax
    formula = ['$log(\textrm{error}) = ', num2str(P(1)), 't+(', num2str(P(2)), ')$'];

    % Create a legend and enable LaTeX interpretation
    legend_entry = legend(formula);
    set(legend_entry, 'Interpreter', 'latex', 'fontsize',13);
    hold off;
end

end
