function [rhohist,psihist, Ham, effSteps,alphathist]=Iter_Fisher(pai,rho0,Q,psi0,alphat,tspan,deltaT,mode)

%%% Initialization
%% Data structure:
% 1) pai,rho0 are row vectors
% 2) y(i,:) is the history of kCur at time i
%    each row of y is a state; each col of y is a time curve of that position

N = length(pai);
TotIt = int64(tspan(2)/deltaT);
rhohist = zeros(TotIt+1,N);            
psihist = zeros(TotIt+1,N);
Ham = zeros(TotIt+1,1);

% effSteps(1) and alphathist(1) are dummy.
effSteps = deltaT*ones(TotIt+1,1);
alphathist =  alphat*ones(TotIt+1,1);

%% Initial Data
rhohist(1,:) = rho0;
psihist(1,:) = psi0;
kCur = rho0./pai;
psiCur = psi0;
% F(p) on table 2 row 4 col 1 
Ham(1) = sum(0.25* pai* (logdiff(kCur).*logdiff(kCur).*Q.*logmean(kCur)));
% first term in eq:Ham-func on Page 4
Ham(1) = Ham(1) + sum(0.25*pai*(logmean(kCur).*psidiffsquare(psiCur).*Q));
% NaN for the dummy value,
effSteps(1) = NaN;          
alphathist(1) = NaN;         
%alpha=-1 means psi is not updated by ODE, the next iteration for rho is MH and the current iteration is not damped.


%%% Iteration

for j=2:(TotIt+1)

    if mod(j,1000)==0
        j
    end

    % If create infeasible value for kCur or psiCur, then end the iteration and save all previous iterations
    if any(kCur<0) || any(isinf(psiCur)) || any(isnan(kCur))
        
        warning('something wrong in (%d)-th iteration',j)
        pause

        rhohist = rhohist(1:(j-1),1:N);
        psihist = psihist(1:(j-1),1:N);
        Ham = Ham(1:(j-1));
        effSteps = effSteps(1:(j-1));
        alphathist = alphathist(1:(j-1));
        break
    end
    
    % warm start by MH, psihist(j)=-log(rhohist(j,:)./pai)
    if deltaT*double(j) < 30
    % if double(j) <= 30
       % if j>2 && Ham(j-1)>Ham(j-2)
       %    warning('Ham increasing')
       %    j
       % end
       rhohist(j,:) = rhohist(j-1,:) + deltaT*(rhohist(j-1,:)*Q);
       kCur = rhohist(j,:)./pai;
       psihist(j,:) = -log(kCur);
       psiCur = psihist(j,:);
       Ham(j) = sum(0.25* pai* (logdiff(kCur).*logdiff(kCur).*Q.*logmean(kCur)));  
       Ham(j) = Ham(j) + sum(0.25*pai*(logmean(kCur).*psidiffsquare(psiCur).*Q));
       alphathist(j) = -1;
       effSteps(j) = deltaT;
       continue
    end


    % design for damping term alpha(t)
    % if deltaT*double(j) > 0
        % alphathist(j) = alphat*log(deltaT*double(j));
    % if deltaT*double(j) >= 30
    % %     %     alphathist(j) = alphat;
    % %     % else
    %         % alphathist(j) = alphat/log(deltaT*double(j));
    % %     % elseif deltaT*double(j) < 500
    % %     %     alphathist(j) = alphat/log(deltaT*double(j));
    % %     % else
    % %     %     alphathist(j) = alphathist(j-1);
    % %     % end
    % %     % alphathist(j) = alphat*log(deltaT*double(j)); 
    %     alphathist(j) = alphat/log(deltaT*double(j));
    %     % if alphathist(j) <= alphat
    %     %     alphathist(j) = alphat;
    %     % end
    % end
    % alphathist(j) = alphat*tanh(deltaT*double(j));
    
    % eq:AMH-kCur on page 8. This equation is the same if one fix logmean.
    LM = logmean(kCur);
    LM = LM.*(-Q);
    LM = RowSumZero(LM);
    LM = LM';
    
    % adjust effective steps
    tmp_deltaT = deltaT;
    negative_part = kCur + psiCur* (tmp_deltaT*LM);
    while any(negative_part(:) < 0)

        % enable in debug mode
        warning('negative part in (%d)-th iteration', j)
        % find(negative_part(:) < 0)
        
        % pause
        tmp_deltaT = 0.1* tmp_deltaT;
        negative_part = kCur + psiCur* (tmp_deltaT*LM);
    end
    effSteps(j) = tmp_deltaT;      
    
    % eq:AMH-kCur on page 8.
    % rho_t = rhohist(j-1, :) + tmp_deltaT * sum(diag(pai) * psidifference_mat(psiCur) .* logmean(rhohist(j-1, :)./pai) .* Q, 2)';
    kCur = kCur + psiCur*(tmp_deltaT*LM);   
    % eq:eQ1 on page 7, replace the log k_i by Table 2 Row 4 Col 2

    %% restart  
    % if it is too far from the target, run MH in the next iteration
    if any(kCur <= 1e-3) %&& deltaT*double(j)< 500
        warning('too far in (%d)-th iteration', j)
        
        psiCur = -log(kCur);
        alphathist(j) = -1;
    else
        % Gauss-Seidel iteration
        psiCur = psiCur + tmp_deltaT*(...
                -alphathist(j)*psiCur...
                - 0.5 * sum((logdiff(kCur)+1-quotient(kCur) + psidiffsquare(psiCur).*partialTheta(kCur)).*Q, 2)'...
               );
        if any(isnan(psiCur))
            j
            find(isnan(psiCur)==1)
        end
    end
    %%
            
    rhohist(j,:) = kCur.*pai;
    psihist(j,:) = psiCur;
    Ham(j) = sum(0.25* pai* (logdiff(kCur).*logdiff(kCur).*Q.*logmean(kCur)));  
    Ham(j) = Ham(j) + sum(0.25*pai*(logmean(kCur).*psidiffsquare(psiCur).*Q));
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
%     title({[num2str(N),'-point convergence in AMH-entropy time stepsize = ', num2str(deltaT, '%0.2e')]
%     ['$alpha_t$=', num2str(alphat, '%0.2e')]
%     });
%     hold off;
% 
    error = sqrt(sum((rhohist(:,1:N)-pai).^2,2));     % cal row norm of the difference
    logError = log10(error);
% 
    subplot(2,1,1)
    % plot(log(t(1:(end-1)/10000:end)),logError(1:(end-1)/10000:end))
    hold on;
    plot(t,logError)
    xlim(tspan)
%     ylim([-3,0])
    % title('M-H: log-log plot of error');
    title({['AMH-entropy: t-log(error) plot']
        ['final timestep=', num2str(deltaT, '%0.2e')]
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
    
    subplot(2,1,2)
    plot(t,Ham)
end


end




% 



