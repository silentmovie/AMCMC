clear;
clc;
close all;

% Get the directory of the current script
scriptDir = fileparts(mfilename('fullpath'));

% Define the relative path to the data file
filename = 'MH-2025-04-03-23-53-52';
parts = split(filename, '-');
dataDir = fullfile(scriptDir, 'data', filename);
parameterFile = fullfile(dataDir, 'parameter.mat');
odeFile = fullfile(dataDir, 'ode.mat');
jumpFile = fullfile(dataDir, 'jump.mat');
paiFile = fullfile(dataDir, 'pai.mat');
% stepFile = fullfile(dataDir, 'steps.mat');
% maxRowFile = fullfile(dataDir, 'maxRowSum.mat');
% One can read the parameter.txt file to setup plot parameter

load(parameterFile)
load(odeFile)
load(jumpFile)
load(paiFile)
% load(stepFile)
% load(maxRowFile)

[TotIt,~] = size(rhoODE);

% TotIt =;

cc = hsv(16);
t = (tspan(1):deltaT:tspan(2))';
    
% pltstep = 0.4/deltaT;         % plt every pltstep iterations
% pltstep = 5*int64(1/deltaT);
pltstep = max(TotIt/100,1);
startpt = 0;
figure('Renderer', 'painters', 'Position', [10 10 1200 1200])
subplot(2,1,1)

L = zeros(3,1);
L(1) = plot(0,nan,'k.');
hold on;
L(2) = plot(0,nan,'k*');
L(3) = plot(0,nan,'k^');

j=0;
for state = 2:4
    j = j+1;
    % plot((startpt:pltstep:TotIt-1),pai(state), 'color', cc(3*(state-1)+1,:), 'marker','.')
    % plot((startpt:pltstep:TotIt-1), rhoODE(startpt+1:pltstep:TotIt, state), 'color', cc(3*(state-1)+1,:),'marker','*');
    % plot((startpt:pltstep:TotIt-1), rhoJump(startpt+1:pltstep:TotIt, state), 'color', cc(3*(state-1)+1,:),'marker','^')
    plot((startpt:pltstep:TotIt-1),pai(state), 'color', cc(2*j-1,:), 'marker','.')
    plot((startpt:pltstep:TotIt-1), rhoODE(startpt+1:pltstep:TotIt, state), 'color', cc(2*j-1,:),'marker','*');
    plot((startpt:pltstep:TotIt-1), rhoJump(startpt+1:pltstep:TotIt, state), 'color', cc(2*j-1,:),'marker','^')
end

hold off;
xlim([startpt,TotIt-1])
% ylim([0,1])
xlabel(['Iteration \times ', num2str(deltaT, '%0.2e'), ' s'])
ylabel('Density')
subtitle({['\fontsize{20} Sampling Dynamics for ', num2str(N), ' nodes']
           ['\fontsize{20} Sampling particles = ', num2str(samplesize, '%0.2e')]
           })
legend(L,'Target',strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',20, 'Location', 'eastoutside');

subplot(2,1,2)

errorODE = sqrt(sum((rhoODE(:,1:N)-pai).^2,2));     % cal row norm of the difference
logErrorODE = log10(errorODE);
errorJump = sqrt(sum((rhoJump(:,1:N)-pai).^2,2));     % cal row norm of the difference
logErrorJump = log10(errorJump);
hold on;
plot((startpt:pltstep:TotIt-1),logErrorODE(startpt+1:pltstep:TotIt),'marker','*')
plot((startpt:pltstep:TotIt-1),logErrorJump(startpt+1:pltstep:TotIt),'marker','^')
plot((startpt:pltstep:TotIt-1),-0.5*log10(samplesize),'Marker','.')
legend(strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',20, 'Location', 'eastoutside');
hold off;
subtitle(['\fontsize{20} t-log(error)'])


% subplot(4,1,2)
% 
% L = zeros(2,1);
% L(1) = plot(0,nan,'k*');
% hold on;
% L(2) = plot(0,nan,'k^');
% 
% for state = 1
%     % plot((0:20), rhoODE(1:21, state), 'color', cc(3*(state-1)+1,:),'marker','*');
%     % plot((0:20), rhoJump(1:21, state), 'color', cc(3*(state-1)+1,:),'marker','^');
%     plot((0:4),pai(state),'Color',cc(state,:),'marker','.');
%     plot((0:4), rhoODE(1:5, state), 'color', cc(state,:),'marker','*');
%     % plot((0:5), rhoJump(1:5, state), 'color', cc(state,:),'marker','^');
% end
% for state = N
%     plot((0:4),pai(state),'Color',cc(state+6-N,:),'marker','.');
%     plot((0:4), rhoODE(1:5, state), 'color', cc(state+6-N,:),'marker','*');
% end
% hold off;
% xlim([0,20])
% % ylim([1/N-0.1,1/N+0.1])
% xlabel('Iteration')
% ylabel('Density')
% subtitle('Validation between ODE and Jump Process for the first few iterations', 'fontsize',15)
% legend(L,strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',13);
% 
% 
% subplot(4,1,2)
% hold on;
% plot((0:pltstep:TotIt-1),HamODE(1:pltstep:TotIt),'marker','*')
% % plot((0:pltstep:TotIt-1),HamJump(1:pltstep:TotIt),'marker','^')
% legend(strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',13,'Location', 'eastoutside');
% hold off;
% subtitle('The Decay of Hamiltonian', 'fontsize',15)


% subplot(4,1,4)
% 
% Diff2MH = sum(log(rhoODE(1:TotIt,:)./pai)+psiODE(1:TotIt,:),2);
% plot((0:TotIt-1),Diff2MH)
% subtitle(['log(p/pi)+psi'])

% subplot(6,1,5)
% hold on;
% plot((0:pltstep:TotIt-2),StepODE(1:pltstep:TotIt-1),'marker','*')
% plot((0:0.5*pltstep:TotIt-2),StepJump(1:0.5*pltstep:TotIt-1),'marker','^')
% legend(strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',13);
% ylim([0.2*deltaT,5*deltaT])
% hold off;
% subtitle('Effective steps', 'fontsize',15)
% 
% subplot(6,1,6)
% hold on;
% plot((0:pltstep:TotIt-2),maxRowODE(1:pltstep:TotIt-1),'marker','*')
% plot((0:0.5*pltstep:TotIt-2),maxRowJump(1:0.5*pltstep:TotIt-1),'marker','^')
% legend(strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',13);
% hold off;
% subtitle('Max Qbar entry in each iteration', 'fontsize',15)

set(findall(gcf, 'Type', 'axes'), 'FontSize', 20);
figFile = fullfile(dataDir, 'Valid-long.png');
saveas(gcf, figFile)