clear;
clc;
close all;

% Get the directory of the current script
scriptDir = fileparts(mfilename('fullpath'));

% Define the relative path to the data file
filename = 'AMHFisher-2025-02-10-11-40-46';
parts = split(filename, '-');
dataDir = fullfile(scriptDir, 'data', filename);
parameterFile = fullfile(dataDir, 'parameter.mat');
odeFile = fullfile(dataDir, 'ode.mat');
jumpFile = fullfile(dataDir, 'jump.mat');
paiFile = fullfile(dataDir, 'pai.mat');
HamFile = fullfile(dataDir, 'ham.mat');
alphatFile = fullfile(dataDir, 'alphat.mat');
stepFile = fullfile(dataDir, 'steps.mat');
% maxRowFile = fullfile(dataDir, 'maxRowSum.mat');
% One can read the parameter.txt file to setup plot parameter

load(parameterFile)
load(odeFile)
load(jumpFile)
load(paiFile)
load(HamFile)
load(alphatFile)
load(stepFile)
% load(maxRowFile)

[TotIt,~] = size(rhoODE);

% TotIt =100;

cc = hsv(6);
t = (tspan(1):deltaT:tspan(2))';
    
% pltstep = 0.4/deltaT;         % plt every pltstep iterations
pltstep = 5*int64(1/deltaT);
% pltstep = 50;
startpt = 0;
figure('Renderer', 'painters', 'Position', [10 10 1200 1200])
%     
subplot(5,1,1)

L = zeros(3,1);
L(1) = plot(0,nan,'k.');
hold on;
L(2) = plot(0,nan,'k*');
L(3) = plot(0,nan,'k^');

j=0;
for state = 1:3
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
subtitle({['\fontsize{18} Sampling Dynamics for ', num2str(N), ' nodes']
           ['\fontsize{18} Sampling particles = ', num2str(samplesize, '%0.2e')]
           ['\fontsize{18} \alpha_t = ', num2str(alphat, '%0.2e')]
           })
legend(L,'Target',strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',18, 'Location', 'eastoutside');


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
subplot(5,1,2)
hold on;
plot((startpt:pltstep:TotIt-1),log10(HamODE(startpt+1:pltstep:TotIt)),'marker','*')
plot((startpt:pltstep:TotIt-1),log10(HamJump(startpt+1:pltstep:TotIt)),'marker','^')
legend(strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',18,'Location', 'eastoutside');
hold off;
subtitle('The Decay of log10(Hamiltonian)', 'fontsize',18)

subplot(5,1,3)

errorODE = sqrt(sum((rhoODE(:,1:N)-pai).^2,2));     % cal row norm of the difference
logErrorODE = log10(errorODE);
errorJump = sqrt(sum((rhoJump(:,1:N)-pai).^2,2));     % cal row norm of the difference
logErrorJump = log10(errorJump);
hold on;
plot((startpt:pltstep:TotIt-1),logErrorODE(startpt+1:pltstep:TotIt),'marker','*')
plot((startpt:pltstep:TotIt-1),logErrorJump(startpt+1:pltstep:TotIt),'marker','^')
plot((startpt:pltstep:TotIt-1),-0.5*log10(samplesize),'Marker','.')
legend(strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',18, 'Location', 'eastoutside');
hold off;
subtitle(['\fontsize{18} t-log(error)'])

subplot(5,1,4)

psisumODE= sum(psiODE,2);
psisumJump = sum(psiJump,2);
hold on;
plot((startpt:pltstep:TotIt-1),psisumODE(startpt+1:pltstep:TotIt),'marker','*')
plot((startpt:pltstep:TotIt-1),psisumJump(startpt+1:pltstep:TotIt),'marker','*')
legend(strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',18, 'Location', 'eastoutside');
hold off;
subtitle(['\fontsize{18} sum of psi'])


% subplot(5,1,5)

% Diff2MHODE = sum(log(rhoODE(1:TotIt,:)./pai)+psiODE(1:TotIt,:),2);
% Diff2MHJump = sum(log(rhoJump(1:TotIt,:)./pai)+psiJump(1:TotIt,:),2);
% hold on;
% plot((startpt:pltstep:TotIt-1),Diff2MHODE(startpt+1:pltstep:TotIt),'marker','*')
% plot((startpt:pltstep:TotIt-1),Diff2MHJump(startpt+1:pltstep:TotIt),'marker','^')
% hold off;
% legend(strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',13, 'Location', 'eastoutside');
% subtitle(['log(p/pi)+psi'])

% subplot(5,1,5)
% plot((0:1200),alphatODE(1:1201))
% subtitle('alphat')

subplot(5,1,5)
hold on;
plot((1:pltstep:TotIt-1),StepODE(2:pltstep:TotIt),'marker','*')
plot((1:0.5*pltstep:TotIt-1),StepJump(2:0.5*pltstep:TotIt),'marker','^')
legend(strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',18);
ylim([0.2*deltaT,5*deltaT])
hold off;
subtitle('Effective steps', 'fontsize',18)
% 
% subplot(6,1,6)
% hold on;
% plot((0:pltstep:TotIt-2),maxRowODE(1:pltstep:TotIt-1),'marker','*')
% plot((0:0.5*pltstep:TotIt-2),maxRowJump(1:0.5*pltstep:TotIt-1),'marker','^')
% legend(strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',13);
% hold off;
% subtitle('Max Qbar entry in each iteration', 'fontsize',15)


figFile = fullfile(dataDir, 'Valid-long.png');
saveas(gcf, figFile)