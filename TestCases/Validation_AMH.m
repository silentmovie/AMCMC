clear;
clc;
close all;

% Get the directory of the current script
scriptDir = fileparts(mfilename('fullpath'));

% Define the relative path to the data file
filename = 'Fisher-2025-05-04-13-40-51';
parts = split(filename, '-');
dataDir = fullfile(scriptDir, 'data', filename);
parameterFile = fullfile(dataDir, 'parameter.mat');
odeFile = fullfile(dataDir, 'ode.mat');
jumpFile = fullfile(dataDir, 'jump.mat');
paiFile = fullfile(dataDir, 'pai.mat');
HamFile = fullfile(dataDir, 'ham.mat');
alphatFile = fullfile(dataDir, 'alphat.mat');
stepFile = fullfile(dataDir, 'steps.mat');

% One can read the parameter.txt file to setup plot parameter

load(parameterFile)
load(odeFile)
load(jumpFile)
load(paiFile)
load(HamFile)
load(alphatFile)
load(stepFile)

pai = pai/(sum(sum(pai)));
[TotIt,~] = size(rhoODE);

% TotIt =650;

cc = hsv(60);
t = (tspan(1):deltaT:tspan(2))';
    
startpt = 0;

% plt every pltstep iterations
pltstep = 1*int64(1/deltaT);
% pltstep = 1;

figure('Renderer', 'painters', 'Position', [10 10 700 1400])
  
% Plot for validation
subplot(3,1,1)

L = zeros(3,1);
L(1) = plot(0,nan,'k.');
hold on;
L(2) = plot(0,nan,'k*');
L(3) = plot(0,nan,'k^');

j=0;
for state = [1,3,64]
    j = j+1;

    plot((startpt:pltstep:TotIt-1),pai(state), 'color', cc(10*j+20,:), 'marker','.')
    plot((startpt:pltstep:TotIt-1), rhoODE(startpt+1:pltstep:TotIt, state), 'color', cc(10*j+20,:),'marker','*');
    plot((startpt:pltstep:TotIt-1), rhoJump(startpt+1:pltstep:TotIt, state), 'color', cc(10*j+20,:),'marker','^')
end

hold off;
% % xlim([startpt,TotIt-1])
% % ylim([0,1])
ylabel('Density')
% subtitle({ ['\fontsize{20} Sampling particles M = ', num2str(samplesize, '%0.2e'), ' distributed in ', num2str(N), ' nodes']
%          })
legend(L,'Target',strcat(parts{1},'-ode'), strcat(parts{1},'-jump'),'fontsize',20, 'Location', 'best');

subplot(3,1,2)
hold on;
plot((startpt:pltstep:TotIt-1),log10(HamODE(startpt+1:pltstep:TotIt)),'b-*')
plot((startpt:pltstep:TotIt-1),log10(HamJump(startpt+1:pltstep:TotIt)),'b-^')
legend(strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',20,'Location', 'best');
hold off;
ylabel('log10(Hamiltonian)')
% subtitle('Decay of log10(Hamiltonian)', 'fontsize',20)

%% plot for error 

subplot(3,1,3)
errorODE = sqrt(sum((rhoODE(:,1:N)-pai).^2,2));     % cal row norm of the difference
logErrorODE = log10(errorODE);
errorJump = sqrt(sum((rhoJump(:,1:N)-pai).^2,2));     % cal row norm of the difference
logErrorJump = log10(errorJump);
hold on;
plot((startpt:pltstep:TotIt-1),logErrorODE(startpt+1:pltstep:TotIt),'marker','*')
plot((startpt:pltstep:TotIt-1),logErrorJump(startpt+1:pltstep:TotIt),'marker','^')
plot((startpt:pltstep:TotIt-1),-0.5*log10(samplesize),'Marker','.')
legend(strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',20, 'Location', 'best');
hold off;
ylabel('log10(error)')
subtitle(['\fontsize{20} t-log(error)'])

%% old plot


% subplot(5,1,4)
% 
% psisumODE= sum(psiODE,2);
% psisumJump = sum(psiJump,2);
% hold on;
% plot((startpt:pltstep:TotIt-1),psisumODE(startpt+1:pltstep:TotIt),'marker','*')
% plot((startpt:pltstep:TotIt-1),psisumJump(startpt+1:pltstep:TotIt),'marker','*')
% legend(strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',20, 'Location', 'eastoutside');
% hold off;
% subtitle(['\fontsize{20} sum of psi'])


% subplot(5,1,5)

% Diff2MHODE = sum(log(rhoODE(1:TotIt,:)./pai)+psiODE(1:TotIt,:),2);
% Diff2MHJump = sum(log(rhoJump(1:TotIt,:)./pai)+psiJump(1:TotIt,:),2);
% hold on;
% plot((startpt:pltstep:TotIt-1),Diff2MHODE(startpt+1:pltstep:TotIt),'marker','*')
% plot((startpt:pltstep:TotIt-1),Diff2MHJump(startpt+1:pltstep:TotIt),'marker','^')
% hold off;
% legend(strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',13, 'Location', 'eastoutside');
% subtitle(['log(p/pi)+psi'])

% subplot(3,1,2)
% hold on;
% plot((1:1:TotIt-1),alphatODE(2:1:TotIt),'b-*')
% plot((1:1:TotIt-1),alphatJump(2:1:TotIt),'b-^')
% legend(strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',20,'Location', 'best');
% hold off;
% subtitle(['\fontsize{20} Damping parameter \gamma(t) = ', num2str(alphat, '%0.2e')])
% 
% subplot(3,1,3)
% hold on;
% plot((1:1:TotIt-1),StepODE(2:1:TotIt),'b-*')
% plot((1:1:TotIt-1),StepJump(2:1:TotIt),'b-^')
% legend(strcat(parts{1},'-','ode'), strcat(parts{1},'-','jump'),'fontsize',20, 'Location', 'best');
% ylim([0.2*deltaT,5*deltaT])
% hold off;
% subtitle('Effective stepsize', 'fontsize',20)
% xlabel(['Iteration \times ', num2str(deltaT, '%0.2e'), ' s'])


set(findall(gcf, 'Type', 'axes'), 'FontSize', 20);

figFile = fullfile(dataDir, 'C3good.png');
saveas(gcf, figFile)