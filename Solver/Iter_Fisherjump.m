function [rhohist, psihist, Ham, effSteps, alphathist] = Iter_Fisherjump(pai, rho0, Q, psi0, alphat, tspan, deltaT, samplesize, mode)

%%% Initialization
%% Data structure:
% 1) pai,rho0 are row vectors
% 2) y(j,:) is the history of k at time j
%    each row of y is a state; each col of y is a time curve of that position

N = length(pai);
TotIt = int64(tspan(2)/deltaT);
rhohist = zeros(TotIt+1, N);
psihist = zeros(TotIt+1,N);
Ham = zeros(TotIt+1, N);
edges = 0.5:1:(N+0.5);           % edges for hiscounts function

% effSteps(1) and alphathist(1) are dummy.
effSteps = deltaT*ones(TotIt+1,1);
alphathist =  0.1141*ones(TotIt+1,1);

%% Initial Data
rhohist(1,:) = rho0;
psihist(1,:) = psi0;
kCur = rho0./pai;
psiCur = psi0;
Ham(1) = sum(0.25* pai* (logdiff(kCur).*logdiff(kCur).*Q.*logmean(kCur)));
Ham(1) = Ham(1) + sum(0.25*pai*(logmean(kCur).*psidiffsquare(psiCur).*Q));
% NaN for the dummy value,
effSteps(1) = NaN;          
alphathist(1) = NaN;         
%alpha=-1 means psi is not updated by ODE, the next iteration for rho is MH and the current iteration is not damped.

% inpo is array of 1*samplesize, 
inpo = randsample(N,samplesize, true, rhohist(1,:))';    % entry is the state of one particle 
% current is array of 1*N
current = histcounts(inpo, N);    % entry is # of particles in each state


%%% Iteration

repeat = 0; % record the number of sampling including repeated resampling

for j = 2:(TotIt+1)

    if mod(j,1000)==0
        j
    end
    
    if any(kCur<0) || any(isinf(psiCur)) || any(isnan(kCur))...
            || any(isnan(psiCur))

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


    % if deltaT*double(j) > 0
    %     % alphathist(j) = alphat*log(deltaT*double(j));
    if deltaT*double(j) >= 30
    %     %     alphathist(j) = alphat;
    %     % else
        % alphathist(j) = alphat/log(deltaT*double(j));
        alphathist(j) =  0.1141/(deltaT*double(j));
        if alphathist(j) <= alphat
            alphathist(j) = alphat;
        end
    %         % alphathist(j) = alphat/log(deltaT*double(j)+30);
    %     % end
    %     % alphathist(j) = alphat*log(deltaT*double(j)); 
    end
    % alphathist(j) = alphat*tanh(deltaT*double(j));

    out = zeros(N,N);
    P = Q.*(PNpsi(psihist(j-1,:)));
    P = diag(pai./rhohist(j-1,:)) * P;
    P = P.*logmean(rhohist(j-1,:) ./ pai );
    P = RowSumZero(P);
    

    tmp_deltaT = deltaT;
    negative_part = eye(N) + tmp_deltaT * P;
    while any(negative_part(:) < 0)
        
        warning('negative part in (%d)-th iteration', j)
        j
        % pause
        tmp_deltaT = 0.1 * tmp_deltaT
        negative_part = eye(N) + tmp_deltaT * P;
    end
    effSteps(j) = tmp_deltaT;
    P = eye(N) + tmp_deltaT * P;
    
    %%% restart method

    % restart by resampling

    tmpcurrent = zeros(1,N);
    while (any (tmpcurrent ==0))
        parfor state = 1:N
            % out is matrix of N*N
            tmp = randsample(N, current(state), true, P(state,:))';
            out(state,:) = histcounts(tmp, edges);
            % out(state,:) = histcounts(tmp, edges.Value);
        end
        tmpcurrent = sum(out,1);
        repeat = repeat + 1;
    end
    repeat = repeat -1;
    current = sum(out,1);
    rhohist(j,:) = current/samplesize;
    kCur = rhohist(j,:)./pai;
    psiCur = psiCur + tmp_deltaT*(...
        -alphathist(j)*psiCur...
        - 0.5 * sum((logdiff(kCur)+1-quotient(kCur)).*Q, 2)'...
        - 0.5* sum(Q.*psidiffsquare(psiCur).*partialTheta(kCur),2)'...
        );

    psihist(j,:) = psiCur;
    Ham(j) = sum(0.25* pai* (logdiff(kCur).*logdiff(kCur).*Q.*logmean(kCur)));  
    Ham(j) = Ham(j) + sum(0.25*pai*(logmean(kCur).*psidiffsquare(psiCur).*Q));
    warning('total repeated sampling = (%d)', repeat)


    % %% restart by MH
    % 
    % parfor state = 1:N
    %     % out is matrix of N*N
    %     tmp = randsample(N, current(state), true, P(state,:))';
    %     out(state,:) = histcounts(tmp, edges);
    %     % out(state,:) = histcounts(tmp, edges.Value);
    % end
    % % futures = cell(N, 1);
    % % for state = 1:N
    % % % futures{state} = parfeval(@computeHist, 1, state, N, current, P, edges);
    % %     futures{state} = parfeval(@computeHist, 1, state, N, current(state), P(state,:), edges);
    % % 
    % % end
    % % 
    % % % Collect results
    % % for state = 1:N
    % %     out(state, :) = fetchOutputs(futures{state});
    % % end
    % 
    % current = sum(out,1);
    % 
    % if any(current==0)
    %     current(current==0) = current(current==0) + 1;
    % 
    %     warning('add particles in (%d)-th iteration, total particle = (%d)', j, sum(current))
    % 
    %     rhohist(j,:) = current/sum(current);
    %     kCur = rhohist(j,:)./pai;
    %     psiCur = -log(kCur);
    %     alphathist(j) = -1;
    % else
    %     rhohist(j,:) = current/sum(current);
    %     kCur = rhohist(j,:)./pai;
    %     % psiCur = psiCur+ tmp_deltaT*(...
    %     %     -alphathist(j)*psiCur...
    %     %     - 0.5 * sum((logdiff(kCur)+1-quotient(kCur)).*Q, 2)'...
    %     %     - 0.5* sum(Q.*psidiffsquare(psiCur).*partialTheta(kCur),2)'...
    %     %     );
    %     psiCur = psiCur + tmp_deltaT*(...
    %     -alphathist(j)*psiCur...
    %     - 0.5 * sum((logdiff(kCur)+1-quotient(kCur) + psidiffsquare(psiCur).*partialTheta(kCur)).*Q, 2)'...
    %             );
    % 
    % 
    % end
    % 
    % psihist(j,:) = psiCur;
    % Ham(j) = sum(0.25* pai* (logdiff(kCur).*logdiff(kCur).*Q.*logmean(kCur)));  
    % Ham(j) = Ham(j) + sum(0.25*pai*(logmean(kCur).*psidiffsquare(psiCur).*Q));

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


    









