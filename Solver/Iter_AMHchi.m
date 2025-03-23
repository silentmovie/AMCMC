function [rhohist,psihist, Ham, effSteps,alphathist] = Iter_AMHchi(pai, rho0, Q, psi0, alphat, tspan, deltaT, mode)

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
    
    cc = hsv(3*N);
    t = [tspan(1):deltaT:tspan(2)]';
    
    pltstep = 20;         % plt every pltstep iterations
    figure('Renderer', 'painters', 'Position', [10 10 1200 1200])
    
    subplot(3,1,1)
    hold on;
    for state = 1:N
        plot((0:TotIt),pai(state), 'color', cc(3*(state-1)+1,:),'marker','.')
        plot((0:pltstep:TotIt), rhohist(1:pltstep:end, state), 'color', cc(3*(state-1)+1,:),'marker','*');
    end

%     xlim([0,tspan(2)])
    ylim([-0.1,1])
    % title('3-point convergence via M-H')
    title({[num2str(N),'-point convergence in AMH-chi time stepsize=', num2str(deltaT, '%0.2e')]
    ['$alpha_t$=', num2str(alphat, '%0.2e')]
    });
    hold off;

    error = sqrt(sum((rhohist(:,1:N)-pai).^2,2));     % cal row norm of the difference
    logError = log10(error);

    subplot(3,1,2)
    % plot(log(t(1:(end-1)/10000:end)),logError(1:(end-1)/10000:end))
    hold on;
    plot(t,error)
    xlim(tspan)
%     ylim([-15,0])
    % title('M-H: log-log plot of error');
    title({['AMH-chi: t-log(error) plot']
        ['timestep=', num2str(deltaT, '%0.2e')]
        });
%     halftime = tspan(2)/2;
%     P = polyfit(t(t>halftime),logError(t>halftime),1);
%     plot(t(t>halftime), P(1)*t(t>halftime) + P(2), 'r', 'LineWidth', 2)
%     % Define the formula using LaTeX syntax
%     formula = ['$log(\textrm{error}) = ', num2str(P(1)), 't+(', num2str(P(2)), ')$'];
% 
%     % Create a legend and enable LaTeX interpretation
%     legend_entry = legend(formula);
%     set(legend_entry, 'Interpreter', 'latex', 'fontsize',13);
%     hold off;
%     
    subplot(3,1,3)
    plot(t,Ham)
    
end


end