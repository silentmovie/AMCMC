clear;
clc;
close all;


%%% Initialization

N = 3;                                        % num of states
seed = 5;                                     % seed for random number generator with method 'twister'
tspan = [0,100];                                % total time span
deltaT = 1e-2;                                % time stepsize
TotIt = tspan(2)/deltaT;                      % total iteration
halftime = tspan(2)/2;
t = [tspan(1):deltaT:tspan(2)]';
mode = 'None';                                 % if 'None', no print, if 'Print', print directly by Iter_AMHchi function

%% create target pai, QMH as Qrow and its minimum eigenvalue
[pai, Qrow, minEig] = ID_Cn(seed, N);

samplesize = 1e6;                              % total particle numbers
% samplesize = ceil(5/min(pai));

%% create initial rho0 and psi0
% rho0 = rand(1,N);
rho0 = ones(1,3);
rho0 = rho0/sum(rho0);

psi0 = -rho0./pai;
% psi0 = rand(1,N)-0.5;
% psi0 = psi0/sum(psi0);
% psi0 = zeros(1,N);



% %% Sample initial data for two-mode guassian
% N = 625;         % num of states

% n = sqrt(N);
% x = linspace(0,1,n);
% y = linspace(0,1,n);
% y = y';

% xc1=0.25;
% yc1=0.25;
% xc2=0.75;
% yc2=0.75;

% % pai = exp(-0.5*((x-xc1).^2+(y-yc1).^2)*20);
% % pai = pai + exp(-0.5*((x-xc2).^2+(y-yc2).^2)*80);
% pai = exp(-0.5*((x-xc1).^2+(y-yc1).^2)*25);

% pai = pai + exp(-0.5*((x-xc2).^2+(y-yc2).^2)*100);
% pai = pai/(sum(sum(pai)));
% pai = reshape(pai,1,N);

% rho0 = ones(1,N);
% rho0 = rho0/sum(rho0);

%% Sample initial data for two-loop

% pai(1,1:40) = 4;
% pai(1,41:60) = 1;
% pai(1,61:100) = 4;
% pai = pai/sum(pai);

% pai=0.2*ones(n,n);
% pai(3,6:19)=1;
% pai(8,6:19)=1;
% pai(17,6:19)=1;
% pai(22,6:19)=1;
% pai(8:17,13)=1;
% pai(3:8,6)=1;
% pai(3:8,19)=1;
% pai(17:22,6)=1;
% pai(17:22,19)=1;

% pai = pai/(sum(sum(pai)));
% pai = reshape(pai,1,N);

% pai = rand(1,N); 
% pai = [2/3, 1/6, 1/6];
% pai = pai/sum(pai);
% pai = sort(pai,'descend');      % target distribution (not necessary to be 'descend')

% rho0 = rand(1,N);
% rho0 = rho0/sum(rho0);

% degree = [2 2 3 2 2 3 2 2];
% degree(1,1:N) = 2;
% degree(1,40) = 3;
% degree(1,61) = 3;           

%% Create a Q-matrix 

% %% Create a Q-matrix and run MH-ode iteration 

% % QMH_{ij}=min{pi_j*deg(i)/[pi_i*deg(j)], 1}
% % QMH(pai, (optional)degree)
% % Q = QMH(pai, degree);
% Q = QMH(pai);
% % edgeQ matrix produces a candidate kernel that embeds graph information,
% % the default one is C_n, each node sends to two nodes uniformly.
% % edge(degree)
% % edge = edgeTwoCycle(degree);
% % edge = edgeCn(pai);
% edge = edgeLattice(pai);
% Q = Q.*edge;              % eq:Q-MH on Page 2

% % RowSumZero is an operation to replace diagonal by each row sum (except on
% % diagonal)
% Qrow = RowSumZero(Q);     % defined between detailed-balanced and Q-MH on Page 2




%% Run MH-method solver 
rhoODE = Iter_MH(pai,rho0,Qrow, tspan,deltaT,mode);

rhoJump = Iter_MHjump(pai, rho0, Qrow, tspan, deltaT, samplesize, mode);

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





