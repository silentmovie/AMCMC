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
    
    cc = hsv(3*N);
    t = [tspan(1):deltaT:tspan(2)]';
    
    pltstep = 1;         % plt every pltstep iterations
    figure('Renderer', 'painters', 'Position', [10 10 1200 1200])
    
%     subplot(2,1,1)
%     hold on;
%     for state = 1:N
%         plot((0:TotIt),pai(state), 'color', cc(3*(state-1)+1,:),'marker','.')
%         plot((0:pltstep:TotIt), rhohist(1:pltstep:end, state), 'color', cc(3*(state-1)+1,:),'marker','*');
%     end
%     
% %     xlim([0,tspan(2)])
%     ylim([-0.1,1])
%     % title('3-point convergence via M-H')
%     title({[num2str(N),'-point convergence time stepsize in MH', num2str(deltaT, '%0.2e')]
%     ['largest negative eigen=', num2str(minEig, '%0.2e')]
%     });
%     hold off;

    error = sqrt(sum((rhohist(:,1:N)-pai).^2,2));     % cal row norm of the difference
    logError = log10(error);

    subplot(2,1,2)
    % plot(log(t(1:(end-1)/10000:end)),logError(1:(end-1)/10000:end))
    hold on;
    plot(t,logError)
    xlim(tspan)
%     ylim([-15,0])
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


% 



