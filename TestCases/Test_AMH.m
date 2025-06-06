clear;
clc;
close all;


%%% Initialization

N = 625;                                        % num of states
seed = 2020;                                     % seed for random number generator with method 'twister'
tspan = [0,1800];                                % total time span
deltaT = 1e-2;                                % time stepsize
TotIt = tspan(2)/deltaT;                      % total iteration
halftime = tspan(2)/2;
t = [tspan(1):deltaT:tspan(2)]';
mode = 'None';                                 % if 'None', no print, if 'Print', print directly by Iter_AMHchi function

method = 'Fisher';                             % method from different theta_{ij} and functional F proposed in the Accelerated method.

%% create target pai, QMH as Qrow and its minimum eigenvalue
% [pai, Qrow, minEigQrow] = ID_Cn(seed, N);
% [pai, Qrow, minEigQrow] = ID_TwoLoop(seed, N);
[pai, Qrow, minEigQrow] =ID_MGaussian(seed, N);
% [pai, Qrow, minEigQrow] =ID_HyperCube(seed, N);

% total particle numbers
% samplesize = ceil(5/min(pai));
samplesize = 5e5;


rayleigh = rayleigh(diag(pai)*Qrow,pai);
% 2*sqrt(-minEigQrow)
% 2*sqrt(rayleigh(1))

%% create initial rho0 and psi0
% rho0 = rand(1,N);
rho0 = ones(1,N);
rho0 = rho0/sum(rho0);

% psi0 = randn(1,N);
psi0 = -log(rho0./pai);
% psi0 = rand(1,N)-0.5;
% psi0 = psi0/sum(psi0);
% psi0 = zeros(1,N);


%% Run AMH method solver 

funODE = str2func(['Iter_',method]);
funJump = str2func(['Iter_',method,'jump']);

for alphat =  2*sqrt(rayleigh(1))     % damping parameter
     
    
    [rhoODE,psiODE,HamODE, StepODE,alphatODE] = funODE(pai,rho0,Qrow,psi0,alphat,tspan,deltaT,mode);
    [rhoJump,psiJump,HamJump, StepJump,alphatJump] = funJump(pai, rho0, Qrow, psi0, alphat, tspan, deltaT, samplesize, seed, mode);
    
    % Auto-save
    Date = datestr(datetime('now'),'yyyy-mm-dd-HH-MM-SS');
    NewFolder = ['data/',method,'-',Date];
    mkdir(NewFolder)

    Qmat = ['Qrow.mat'];
    Qmat = fullfile(NewFolder, Qmat);

    ODEmat = ['ode.mat'];
    ODEmat = fullfile(NewFolder, ODEmat);

    Jumpmat = ['jump.mat'];
    Jumpmat = fullfile(NewFolder, Jumpmat);

    paimat = ['pai.mat'];
    paimat = fullfile(NewFolder, paimat);

    Hammat = ['ham.mat'];
    Hammat = fullfile(NewFolder, Hammat);

    parametermat =  ['parameter.mat'];
    parametermat = fullfile(NewFolder, parametermat);
   
    
    stepmat = ['steps.mat'];
    stepmat = fullfile(NewFolder, stepmat);

    alphatmat = ['alphat.mat'];
    alphatmat = fullfile(NewFolder, alphatmat);

    save(Qmat, 'Qrow');
    save(ODEmat, 'rhoODE', 'psiODE');
    save(paimat, 'pai');
    save(parametermat, 'N', 'seed', 'tspan', 'deltaT', 'samplesize','alphat','minEigQrow');
    
    save(Jumpmat, 'rhoJump', 'psiJump');
    save(Hammat, 'HamODE', 'HamJump');
    save(alphatmat,'alphatODE','alphatJump');
    save(stepmat, 'StepODE', 'StepJump');


    parametertxt = ['parameter.txt'];
    parametertxt = fullfile(NewFolder, parametertxt);

    parameterlist = {
        ['N = ' num2str(N)]
        ['seed = ' num2str(seed)]
        ['tspan = ' num2str(tspan)]
        ['deltaT = ' num2str(deltaT, '%0.2e')]
        ['samplesize = ' num2str(samplesize, '%0.2e')]
        ['end_alphat = ' num2str(alphat, '%0.2e')]
        ['minEigQrow = ' num2str(minEigQrow, '%0.2e')]
        };

    fid = fopen(parametertxt,'w');
    fmtString = [repmat('%s\t',1,size(parameterlist,2)-1),'%s\n'];
    fprintf(fid,fmtString,parameterlist{:});
    fclose(fid);

end




