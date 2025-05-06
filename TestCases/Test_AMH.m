clear;
clc;
close all;


%%% Initialization

N = 64;                                        % num of states
seed = 2596;                                     % seed for random number generator with method 'twister'
tspan = [0,100];                                % total time span
deltaT = 1e-2;                                % time stepsize
TotIt = tspan(2)/deltaT;                      % total iteration
halftime = tspan(2)/2;
t = [tspan(1):deltaT:tspan(2)]';
mode = 'None';                                 % if 'None', no print, if 'Print', print directly by Iter_AMHchi function

method = 'Fisher';                             % method from different theta_{ij} and functional F proposed in the Accelerated method.

%% create target pai, QMH as Qrow and its minimum eigenvalue
% [pai, Qrow, minEig] = ID_Cn(seed, N);
% [pai, Qrow, minEig] = ID_TwoLoop(seed, N);
% [pai, Qrow, minEigQrow] =ID_MGaussian(seed, N);
[pai, Qrow, minEigQrow] =ID_HyperCube(seed, N);
% samplesize = 484314;                              % total particle numbers
% samplesize = ceil(5/min(pai));
samplesize = 1e4;
% 2*sqrt(-minEigQrow)
% rayleigh = rayleigh(-diag(pai)*Qrow,pai);
rayleigh = rayleigh(diag(pai)*Qrow,pai);
2*sqrt(rayleigh(1))
pause
%% create initial rho0 and psi0
% rho0 = rand(1,N);
rho0 = ones(1,N);
rho0 = rho0/sum(rho0);
% rho0 =[0.138558909613056	0.138558909613056	0.150907613766526	0.0719745670073633	0.0719745670073633	0.150907613766526	0.138558909613056	0.138558909613056];
% psi0 =[0.0669171987353261	0.0669171987353261	-0.0184550461125071	-0.258929299518705	-0.258929299518705	-0.0184550461125071	0.0669171987353260	0.0669171987353260];
% psi0 = randn(1,N);
psi0 = -log(rho0./pai);
% psi0 = rand(1,N)-0.5;
% psi0 = psi0/sum(psi0);
% psi0 = zeros(1,N);


%% Run AMH method solver 

funODE = str2func(['Iter_',method]);
funJump = str2func(['Iter_',method,'jump']);
2*sqrt(rayleigh(1))
for alphat = 0.17       % damping parameter
     
    
    [rhoODE,psiODE,HamODE, StepODE,alphatODE] = funODE(pai,rho0,Qrow,psi0,alphat,tspan,deltaT,mode);
    [rhoJump,psiJump,HamJump, StepJump,alphatJump] = funJump(pai, rho0, Qrow, psi0, alphat, tspan, deltaT, samplesize, mode);
    
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
    % 
    % 
    % save(Hammat, 'HamODE')
    % save(alphatmat,'alphatODE')
    % save(stepmat, 'StepODE');

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




