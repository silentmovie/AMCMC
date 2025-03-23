% clear;
% clc;
% close all;


%% Initialization
rng(4,"twister")
N = 3;         % num of states
% pai = rand(1,N);                   
% pai = pai/sum(pai);
% pai = sort(pai,'descend');      % target distribution (not necessary to be 'descend'
% pai = ones(1,N);
% pai = pai/sum(pai);

pai = [2/3, 1/6, 1/6];
rho0 = rand(1,N);
rho0 = rho0/sum(rho0);

% rho0 = [0.1, 0.1, 0.8]

% psi0 = rand(1,N)-0.5;
% psi0 = psi0/sum(psi0);
% psi0 = zeros(1,N);
psi0 = -rho0./pai;

% alphat = sqrt(3);

tspan = [0,50];
deltaT = 1e-1;
TotIt = tspan(2)/deltaT;
halftime = tspan(2)/2;
t = [tspan(1):deltaT:tspan(2)]';
mode = 'None';


%% Create a Q-matrix and run MH-ode iteration 

Q = QMH(pai);
edge = edgeQ(pai);
Q = Q.*edge;
Qrow = RowSumZero(Q)

Eigenvalue = eig(Qrow);
% minEig = max(Eigenvalue(abs(Eigenvalue)>=1e-3));
% alphat = 2*sqrt(-minEig);
% alphat = 0;

%% Run AMH-chi method solver 

for alphat = 2*sqrt(-minEig)        % damping parameter
     

    [rhoODE,psiODE,HamODE, StepODE,alphatODE] = Iter_AMHchi(pai,rho0,Qrow,psi0,alphat,tspan,deltaT,mode);
    [rhoJump,psiJump,HamJump, StepJump,alphatJump] = Iter_AMHchijump(pai, rho0, Qrow, psi0, alphat, tspan, deltaT, samplesize, mode);

    % Auto-save
    Date = datestr(datetime('now'),'yyyy-mm-dd-HH-MM-SS');
    NewFolder = ['data/AMHFisher-',Date];
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
    save(parametermat, 'N', 'seed', 'tspan', 'deltaT', 'samplesize','alphat','minEig');
    
    % save(Jumpmat, 'rhoJump', 'psiJump');
    % save(Hammat, 'HamODE', 'HamJump');
    % save(alphatmat,'alphatODE','alphatJump');
    % save(stepmat, 'StepODE', 'StepJump');
    
    % 
    save(Hammat, 'HamODE')
    save(alphatmat,'alphatODE')
    save(stepmat, 'StepODE');

    parametertxt = ['parameter.txt'];
    parametertxt = fullfile(NewFolder, parametertxt);

    parameterlist = {
        ['N = ' num2str(N)]
        ['seed = ' num2str(seed)]
        ['tspan = ' num2str(tspan)]
        ['deltaT = ' num2str(deltaT, '%0.2e')]
        ['samplesize = ' num2str(samplesize, '%0.2e')]
        ['ini_alphat = ' num2str(alphat, '%0.2e')]
        ['minEig = ' num2str(minEig, '%0.2e')]
        };

    fid = fopen(parametertxt,'w');
    fmtString = [repmat('%s\t',1,size(parameterlist,2)-1),'%s\n'];
    fprintf(fid,fmtString,parameterlist{:});
    fclose(fid);

end




