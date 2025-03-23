function [rhohist, psihist, Ham, effSteps, alphathist] = Iter_AMHFisherjump(pai, rho0, Q, psi0, alphat, tspan, deltaT, samplesize, mode)

N = length(pai);

TotIt = int64(tspan(2)/deltaT);
rhohist = zeros(TotIt+1, N);
psihist = zeros(TotIt+1, N);
Ham = zeros(TotIt+1, 1);
effSteps = zeros(TotIt+1,1);
alphathist = alphat*ones(TotIt+1,1);

edges = 0.5:1:(N+0.5);


% PMH = eye(N) + Q*deltaT;          %eq: MH-jump on Page 3

% inpo is array of 1*samplesize, 
inpo = randsample(N,samplesize, true, rho0)';    % entry is the state of one particle 
% current is array of 1*N
current = histcounts(inpo, N);    % entry is # of particles in each state


rhohist(1,:) = rho0;
psihist(1,:) = psi0;
psi = psi0;
k = rho0./pai;

Ham(1) = sum(0.25* pai* (logdiff(k).*logdiff(k).*Q.*logmean(k)));
Ham(1) = Ham(1) + sum(0.25*pai*(logmean(k).*psidiffsquare(psi0).*Q));

% tic
for j = 2:(TotIt+1)

    if mod(j,1000)==0
        j
        % toc 
        % tic
    end
    
    if any(k<0) || any(isinf(psi)) || any(isnan(k))...
            || any(isnan(psi))

        warning('something wrong')
        j
        pause

        rhohist = rhohist(1:(j-1),1:N);
        psihist = psihist(1:(j-1),1:N);
        Ham = Ham(1:(j-1));
        effSteps = effSteps(1:(j-1));
        alphathist = alphathist(1:(j-1));
        break
    end


    out = zeros(N,N);
    
    if deltaT*double(j) <= 0
    % if double(j) <= 0
        
       rhohist(j,:) = rhohist(j-1,:) + deltaT*(rhohist(j-1,:)*Q);
       k = rhohist(j,:)./pai;
       psihist(j,:) = -log(k);
       psi = psihist(j,:);
       Ham(j) = sum(0.25* pai* (logdiff(k).*logdiff(k).*Q.*logmean(k)));  
       Ham(j) = Ham(j) + sum(0.25*pai*(logmean(k).*psidiffsquare(psi).*Q));
       alphathist(j) = -1;
       effSteps(j) = deltaT;
       continue
    end


    % if Ham(j-1) > Ham(j-2)
    % 
    %    warning('Ham increase')
    %    j
    %    parfor state = 1:N
    %         %P(state,:) is the transition prob from state to the other probability
    % 
    %         tmp = randsample(N, current(state), true, PMH(state,:))';
    %         out(state,:) = histcounts(tmp, edges);
    %    end
    %    current = sum(out,1);    % sum along column
    % 
    %    rhohist(j,:) = current/samplesize;
    %    k = rhohist(j,:)./pai;
    %    psihist(j,:) = -log(k);
    %    psi = psihist(j,:);
    %    Ham(j) = sum(0.25* pai* (logdiff(k).*logdiff(k).*Q.*logmean(k)));
    %    Ham(j) = Ham(j) + sum(0.25*pai*(logmean(k).*psidiffsquare(psi).*Q));
    %    alphathist(j) = 0;
    %    continue
    % end
    if deltaT*double(j) > 0
        % alphathist(j) = alphat*log(deltaT*double(j));
        % if deltaT*double(j) < exp(1)
        %     alphathist(j) = alphat;
        % else
            alphathist(j) = alphat/log(deltaT*double(j)+30);
        % end
        % alphathist(j) = alphat*log(deltaT*double(j)); 
    end
    % alphathist(j) = alphat;
    
    P = Q.*(PNpsi(psihist(j-1,:)));
    P = diag(pai./rhohist(j-1,:)) * P;
    P = P.*logmean(rhohist(j-1,:) ./ pai );
    P = RowSumZero(P);
    

    tmp_deltaT = deltaT;
    negative_part = eye(N) + tmp_deltaT * P;
    while any(negative_part(:) < 0)
        
        warning('negative part')
        j
        % pause
        tmp_deltaT = 0.1 * tmp_deltaT
        negative_part = eye(N) + tmp_deltaT * P;
    end
    effSteps(j) = tmp_deltaT;
    P = eye(N) + tmp_deltaT * P;

    % if 4/(deltaT*N) < min(min(abs(P)))
    %     min(min(abs(P)))
    %     pause
    % end

    
    parfor state = 1:N
        % out is matrix of N*N
        tmp = randsample(N, current(state), true, P(state,:))';
        out(state,:) = histcounts(tmp, edges);
        % out(state,:) = histcounts(tmp, edges.Value);
    end
    
    % current(i) = # of particle at node i.
    current = sum(out,1);
    
    %% !!Artificial Approximation together
    if any(current==0)
        current(current==0) = current(current==0) + 1;
        
        warning('add particles')
        j
        if sum(current) ~= samplesize
            warning(['current particle = ',num2str(sum(current))])
        end
        % -alphat*psi
        % -0.5 * sum((logdiff(k)+1-quotient(k)).*Q, 2)'
        % -0.5 * sum(Q.*psidiffsquare(psi).*omega(k),2)'
        % pause
       
	    rhohist(j,:) = current/sum(current);
        % current = floor(rhohist(j,:)*samplesize);
        k = rhohist(j,:)./pai;
        psi = -log(k);
        alphathist(j) = -1;
    else
        rhohist(j,:) = current/sum(current);
        k = rhohist(j,:)./pai;
        psi = psi + tmp_deltaT*(...
            -alphathist(j)*psi...
            - 0.5 * sum((logdiff(k)+1-quotient(k)).*Q, 2)'...
            - 0.5* sum(Q.*psidiffsquare(psi).*omega(k),2)'...
            );


    end
    %%
    
    
    psihist(j,:) = psi;
    Ham(j) = sum(0.25* pai* (logdiff(k).*logdiff(k).*Q.*logmean(k)));  
    Ham(j) = Ham(j) + sum(0.25*pai*(logmean(k).*psidiffsquare(psi).*Q));
    
%     psi = psi + tmp_deltaT*(...
%             -alphat*psi...
%             -sum((logdiff(k)+1-quotient(k)).*Q, 2)'...
%             - 0.5* sum(Q.*psidiffsquare(psi).*omega(k),2)'...
%             );

    % badterm = (1- quotient(k)).* (1 - psidiffsquare(psi)./(logdiff(k).*logdiff(k)));
    % psi = psi + tmp_deltaT*(...
    %         -alphat*psi...
    %         - 0.5 * sum(Q.* (logdiff(k) + psidiffsquare(psi)./logdiff(k) + badterm) ,2)'...
    %         );
    %fprintf(" Value of psi eq right: ");
  

end
% toc

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
% %     title({[num2str(N),'-point convergence time stepsize ', num2str(deltaT, '%0.2e')]
% %     ['largest negative eigen=', num2str(minEig, '%0.2e')]
% %     });
%     hold off;

    error = sqrt(sum((rhohist(:,1:N)-pai).^2,2));     % cal row norm of the difference
    logError = log10(error);

    subplot(2,1,2)
    % plot(log(t(1:(end-1)/10000:end)),logError(1:(end-1)/10000:end))
    hold on;
    plot(t,logError)
%     xlim(tspan)
%     ylim([-15,0])
%     % title('M-H: log-log plot of error');
%     title({['MH: t-log(error) plot']
%         ['timestep=', num2str(deltaT, '%0.2e')]
%         });
    halftime = tspan(2)/2;
    P = polyfit(t(t>halftime),logError(t>halftime),1);
    plot(t(t>halftime), P(1)*t(t>halftime) + P(2), 'r', 'LineWidth', 2)
    % Define the formula using LaTeX syntax
    formula = ['$log(\textrm{error}) = ', num2str(P(1)), 't+(', num2str(P(2)), ')$'];
% 
    % Create a legend and enable LaTeX interpretation
    legend_entry = legend(formula);
    set(legend_entry, 'Interpreter', 'latex', 'fontsize',13);
    hold off;
end

end


    









