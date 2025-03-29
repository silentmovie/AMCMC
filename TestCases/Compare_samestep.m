clear;
clc;
close all;

% Get the directory of the current script
scriptDir = fileparts(mfilename('fullpath'));

MHfolder = 'MH-2025-03-28-23-46-35-twocycle-figure3';
AMHfolder = 'Fisher-2025-03-29-00-39-32';

dataDir = fullfile(scriptDir, 'data', MHfolder);
parameterFile = fullfile(dataDir, 'parameter.mat');
MHodeFile = fullfile(dataDir, 'ode.mat');
MHjumpFile = fullfile(dataDir, 'jump.mat');
paiFile = fullfile(dataDir, 'pai.mat');

load(parameterFile)
load(MHodeFile)
load(MHjumpFile)
load(paiFile)

% TotIter is the total iterations generated from Iter_solver,
% TotIt is the total iterations used for plot

[TotIter,~] = size(rhoODE);
TotIt = TotIter;
startpt = 0;
pltstep = max(TotIt/200,1);
% pltstep = 0.4/deltaT;         % plt every pltstep iterations
% pltstep = 5*int64(1/deltaT);



figure('Renderer', 'painters', 'Position', [10 10 600 400])
hold on;

errorODE = sqrt(sum((rhoODE(:,1:N)-pai).^2,2));     % cal row norm of the difference
logErrorODE = log10(errorODE);
errorJump = sqrt(sum((rhoJump(:,1:N)-pai).^2,2));     % cal row norm of the difference
logErrorJump = log10(errorJump);

% xlim([0,100])
xlabel(['\fontsize{20} Iteration \times ', num2str(deltaT, '%0.2e'), ' s'])
plot((0:pltstep:TotIt-1),logErrorODE(1:pltstep:TotIt),'r-*')
plot((0:pltstep:TotIt-1),logErrorJump(1:pltstep:TotIt),'r-^')
% LegendPlot{1}=legend({'$\frac{1}{\sqrt{\text{sample size}}}$'}, 'Interpreter', 'latex');
% LegendPlot = cell(5,1);
% LegendPlot{1} = '\frac{1}{\sqrt{\text{sample size}}}';
% LegendPlot{2} = 'MH-ode';
% LegendPlot{3} = 'MH-jump';


dataDir = fullfile(scriptDir, 'data', AMHfolder);
parts = split(AMHfolder, '-');
parameterFile = fullfile(dataDir, 'parameter.mat');
AMHodeFile = fullfile(dataDir, 'ode.mat');
AMHjumpFile = fullfile(dataDir, 'jump.mat');
HamFile = fullfile(dataDir, 'ham.mat');
alphatFile = fullfile(dataDir, 'alphat.mat');
stepsFile = fullfile(dataDir, 'steps.mat');

load(parameterFile)
load(AMHodeFile)
load(AMHjumpFile)
load(HamFile)
load(alphatFile)
load(stepsFile)


errorODE = sqrt(sum((rhoODE(:,1:N)-pai).^2,2));     % cal row norm of the difference
logErrorODE = log10(errorODE);
errorJump = sqrt(sum((rhoJump(:,1:N)-pai).^2,2));     % cal row norm of the difference
logErrorJump = log10(errorJump);


plot((0:pltstep:TotIt-1),logErrorODE(1:pltstep:TotIt),'b-*')
plot((0:pltstep:TotIt-1),logErrorJump(1:pltstep:TotIt),'b-^')

% LegendPlot = cell(5,1);
% LegendPlot{1} = '\frac{1}{\sqrt{\text{sample size}}}';
% LegendPlot{2} = 'MH-ode';
% LegendPlot{3} = 'MH-jump';
% LegendPlot{4} = [parts{1},'-ode'];
% LegendPlot{5} = [parts{1},'-jump'];

plot((0:pltstep:TotIt-1),-0.5*log10(samplesize) + 0.0*(0:pltstep:TotIt-1),'k-')
plot((0:pltstep:TotIt-1),-log10(samplesize) + 0.0*(0:pltstep:TotIt-1),'k-.')


LegendPlot = legend('\textrm{MH-ode}','\textrm{MH-jump}',  ...
                    [parts{1},'-ode'],[parts{1},'-jump'], ...
                    '$\log\frac{1}{\sqrt{M}}$','$\log\frac{1}{M}$');
set(LegendPlot,'Interpreter','latex','fontsize',20, 'Location', 'best'); %'eastoutside'
set(LegendPlot,'FontName', 'Helvetica')


% text(-3,-0.5*log10(samplesize),,'Interpreter', 'latex', 'FontSize', 20);
% text(-1.5,-log10(samplesize),'$\log\frac{1}{\sqrt{M}}$','Interpreter', 'latex', 'FontSize', 20);
hold off;


set(findall(gcf, 'Type', 'axes'), 'FontSize', 24);
figFile = fullfile(dataDir, 'Compare.png');
saveas(gcf, figFile)



