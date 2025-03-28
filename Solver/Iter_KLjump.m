function [rhohist, Ham] = Iter_KLjump(pai, rho0, Q, psi0, alphat, tspan, deltaT, samplesize, mode)

N = length(pai);

TotIt = int64(tspan(2)/deltaT);
rhohist = zeros(TotIt+1, N);
Ham = zeros(TotIt+1,1);
% inpo is array of 1*samplesize, 
inpo = randsample(N,samplesize, true, rho0)';    % entry is the state of one particle 
% current is array of 1*N
current = histcounts(inpo, N);    % entry is # of particles in each state
% out is matrix of N*N
edges = [0.5:1:(N+0.5)];

rhohist(1,:) = rho0;
psi = psi0;
k = rho0./pai;

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
        tmp = randsample(N, current(state), true, W(state,:))';
        out(state,:) = histcounts(tmp, edges);
    end
    
    current = sum(out,1);
    rhohist(j,:) = current/samplesize;
    psi = psi + deltaT*(-alphat*psi -log(k) - 0.5* sum(Q.*psidiffsquare(psi).*omega(k),2)');
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
    
    cc = hsv(3*N);
    t = [tspan(1):deltaT:tspan(2)]';
    
    pltstep = 1;         % plt every pltstep iterations
    figure('Renderer', 'painters', 'Position', [10 10 1200 1200])
    
    subplot(2,1,1)
    hold on;
    for state = 1:N
        plot((0:TotIt),pai(state), 'color', cc(3*(state-1)+1,:),'marker','.')
        plot((0:pltstep:TotIt), rhohist(1:pltstep:end, state), 'color', cc(3*(state-1)+1,:),'marker','*');
    end
    
%     xlim([0,tspan(2)])
    ylim([-0.1,1])
    % title('3-point convergence via M-H')
%     title({[num2str(N),'-point convergence time stepsize ', num2str(deltaT, '%0.2e')]
%     ['largest negative eigen=', num2str(minEig, '%0.2e')]
%     });
    hold off;

    error = sqrt(sum((rhohist(:,1:N)-pai).^2,2));     % cal row norm of the difference
    logError = log10(error);

%     subplot(2,1,2)
%     % plot(log(t(1:(end-1)/10000:end)),logError(1:(end-1)/10000:end))
%     hold on;
%     plot(t,logError)
%     xlim(tspan)
%     ylim([-15,0])
%     % title('M-H: log-log plot of error');
%     title({['MH: t-log(error) plot']
%         ['timestep=', num2str(deltaT, '%0.2e')]
%         });
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
end

end

    