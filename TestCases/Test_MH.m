clear;
clc;
close all;


%%% Initialization

N = 625;                                        % num of states
seed = 2596;                                     % seed for random number generator with method 'twister'
tspan = [0,3000];                                % total time span
deltaT = 1e-2;                                % time stepsize
TotIt = tspan(2)/deltaT;                      % total iteration
halftime = tspan(2)/2;
t = [tspan(1):deltaT:tspan(2)]';
mode = 'None';                                 % if 'None', no print, if 'Print', print directly by Iter_AMHchi function

%% create target pai, QMH as Qrow and its minimum eigenvalue
% [pai, Qrow, minEig] = ID_Cn(seed, N);
% [pai, Qrow, minEig] = ID_TwoCycle(seed, N);
[pai, Qrow, minEig] = ID_MGaussian(seed, N);
% [pai, Qrow, minEigQrow] = ID_HyperCube(seed, N);


% total particle numbers
% samplesize = ceil(5/min(pai));
samplesize = 5e5;
%% create initial rho0
% rho0 = rand(1,N);
rho0 = ones(1,N);
rho0 = rho0/sum(rho0);


%% Run MH-method solver 
rhoODE = Iter_MH(pai,rho0,Qrow, tspan,deltaT,mode);

rhoJump = Iter_MHjump(pai, rho0, Qrow, tspan, deltaT, samplesize, seed,mode);

%% Auto-save
Date = datestr(datetime('now'),'yyyy-mm-dd-HH-MM-SS');
NewFolder = ['data/MH-',Date];
mkdir(NewFolder)

ODEmat = ['ode.mat'];
ODEmat = fullfile(NewFolder, ODEmat);

Jumpmat = ['jump.mat'];
Jumpmat = fullfile(NewFolder, Jumpmat);

paimat = ['pai.mat'];
paimat = fullfile(NewFolder, paimat);

parametermat =  ['parameter.mat'];
parametermat = fullfile(NewFolder, parametermat);


save(ODEmat, 'rhoODE');
save(Jumpmat, 'rhoJump');
save(paimat, 'pai');
save(parametermat, 'N', 'seed', 'tspan', 'deltaT', 'samplesize');


parametertxt = ['parameter.txt'];
parametertxt = fullfile(NewFolder, parametertxt);

parameterlist = {
    ['N = ' num2str(N)]
    ['seed = ' num2str(seed)]
    ['tspan = ' num2str(tspan)]
    ['deltaT = ' num2str(deltaT, '%0.2e')]
    ['samplesize = ' num2str(samplesize, '%0.2e')]
    ['minEig = ' num2str(minEig, '%0.2e')]
    };

fid = fopen(parametertxt,'w');
fmtString = [repmat('%s\t',1,size(parameterlist,2)-1),'%s\n'];
fprintf(fid,fmtString,parameterlist{:});
fclose(fid);





