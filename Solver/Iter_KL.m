function [rhohist,psihist, Ham, effSteps,alphathist]=Iter_KL(pai,rho0,Q,psi0,alphat,tspan,deltaT,mode)

%%% Initialization
%% Data structure:
% 1) pai,rho0 are row vectors
% 2) y(j,:) is the history of k at time j
%    each row of y is a state; each col of y is a time curve of that position

N = length(pai);
TotIt = int64(tspan(2)/deltaT);
rhohist = zeros(TotIt+1,N);            %each row is a state; each col is a time curve of that position
psihist = zeros(TotIt+1,N);
Ham = zeros(TotIt+1,1);

% effSteps(1) and alphathist(1) are dummy.
effSteps = deltaT*ones(TotIt+1,1);
alphathist = alphat*ones(TotIt+1,1);

%% Initial Data
rhohist(1,:) = rho0;               
psihist(1,:) = psi0;
kCur = rho0./pai;
psiCur = psi0;

Ham(1) = sum(log(rhohist(1,:)).*rhohist(1,:) - log(pai).*rhohist(1,:));
Ham(1) = Ham(1) + sum(0.25*pai*(logmean(kCur).*psidiffsquare(psi0).*Q));

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

    % warm-stary by MH, psihist(j)=-log(rhohist(j,:)./pai)
    if deltaT*double(j) <=0
    % if double(j) <= 30
       % if j>2 && Ham(j-1)>Ham(j-2)
       %    warning('Ham increasing')
       %    j
       % end
        rhohist(j,:) = rhohist(j-1,:) + deltaT*(rhohist(j-1,:)*Q);
        kCur = rhohist(j,:)./pai;
        psihist(j,:) = -log(kCur);
        psiCur = psihist(j,:);
        Ham(j) = sum(log(rhohist(j,:)).*rhohist(j,:) - log(pai).*rhohist(j,:));
        Ham(j) = Ham(j) + sum(0.25*pai*(logmean(kCur).*psidiffsquare(psiCur).*Q));
        alphathist(j) = -1;
        effSteps(j) = deltaT;
        continue
    end


    % design for damping term alpha(t)
    if deltaT*double(j) > 0
        % alphathist(j) = alphat*log(deltaT*double(j));
        % if deltaT*double(j) < exp(1)
        %     alphathist(j) = alphat;
        % else
            alphathist(j) = alphat/log(deltaT*double(j)+30);
        % elseif deltaT*double(j) < 500
        %     alphathist(j) = alphat/log(deltaT*double(j));
        % else
        %     alphathist(j) = alphathist(j-1);
        % end
        % alphathist(j) = alphat*log(deltaT*double(j)); 
    end
    % alphathist(j) = alphat;
    

    LM = logmean(kCur);
    LM = LM.*(-Q);
    LM = RowSumZero(LM);
    LM = LM';


    % adjust effective steps
    tmp_deltaT = deltaT;
    negative_part = kCur + psiCur* (tmp_deltaT*LM);
    while any(negative_part(:) < 0)

        % enable in debug mode
        % warning('negative part in (%d)-th iteration', j)
        % find(negative_part(:) < 0)
        
        % pause
        tmp_deltaT = 0.1* tmp_deltaT;
        negative_part = kCur + psiCur* (tmp_deltaT*LM);
    end
    effSteps(j) = tmp_deltaT;      
    
    % eq:AMH-kCur on page 8.
    % rho_t = rhohist(j-1, :) + tmp_deltaT * sum(diag(pai) * psidifference_mat(psiCur) .* logmean(rhohist(j-1, :)./pai) .* Q, 2)';
    kCur = kCur + psiCur*(tmp_deltaT*LM);  
  


    %% restart
    % if it is too far from the target, run MH in the next iteration
    if any(kCur <= 1e-3) %&& deltaT*double(j)< 500
        warning('too far in (%d)-th iteration', j)
        
        psiCur = -log(kCur);
        alphathist(j) = -1;
    else
        % Gauss-Seidel iteration
        psiCur = psiCur + deltaT*(-alphathist(j)*psiCur -log(kCur) - 0.5* sum(Q.*psidiffsquare(psiCur).*partialTheta(kCur),2)');
%     psi = psi + deltaT*(-alphat*psi -2*pai*logdiff(tmp)*Q./tmp - 0.5* sum(Q.*psidiffsquare(psi).*partialTheta(tmp),2)');
        if any(isnan(psiCur))
            j
            find(isnan(psiCur)==1)
        end
    end

    rhohist(j,:) = kCur.*pai;
    psihist(j,:) = psiCur;
    Ham(j) = sum(log(rhohist(j,:)).*rhohist(j,:) - log(pai).*rhohist(j,:));
    Ham(j) = Ham(j) + sum(0.25*pai*(logmean(kCur).*psidiffsquare(psiCur).*Q));
end


if strcmp(mode, 'None')
    
    return
    
elseif strcmp(mode, 'Print')
    
   warning('use plot function in TestCases folder')
   return
end


end




% 



