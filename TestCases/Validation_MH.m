clear;
clc;
close all;

% Get the directory of the current script
scriptDir = fileparts(mfilename('fullpath'));

% Define the relative path to the data file
filename = 'MH-2025-05-10-10-18-20';
parts = split(filename, '-');
dataDir = fullfile(scriptDir, 'data', filename);
parameterFile = fullfile(dataDir, 'parameter.mat');
odeFile = fullfile(dataDir, 'ode.mat');
jumpFile = fullfile(dataDir, 'jump.mat');
paiFile = fullfile(dataDir, 'pai.mat');
% stepFile = fullfile(dataDir, 'steps.mat');
% One can read the parameter.txt file to setup plot parameter

load(parameterFile)
load(odeFile)
load(jumpFile)
load(paiFile)
% load(stepFile)

pai = pai/(sum(sum(pai)));
[TotIt,~] = size(rhoODE);


cc = hsv(60);
t = (tspan(1):deltaT:tspan(2))';
    
% plt every pltstep iterations
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
for state = [1,3,64]
    j = j+1;
    plot((startpt:pltstep:TotIt-1),pai(state), 'color', cc(10*j+20,:), 'marker','.')
    plot((startpt:pltstep:TotIt-1), rhoODE(startpt+1:pltstep:TotIt, state), 'color', cc(10*j+20,:),'marker','*');
    plot((startpt:pltstep:TotIt-1), rhoJump(startpt+1:pltstep:TotIt, state), 'color', cc(10*j+20,:),'marker','^')
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

%% old plot


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


set(findall(gcf, 'Type', 'axes'), 'FontSize', 20);
figFile = fullfile(dataDir, 'Valid-long.png');
saveas(gcf, figFile)