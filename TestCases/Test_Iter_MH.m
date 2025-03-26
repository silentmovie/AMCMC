% clear;
% clc;
% close all;


%% Initialization
rng(8,"twister")
N = 10;  % 5;  % 10;         % num of states
pai = rand(1,N);                   
pai = pai/sum(pai);
pai = sort(pai,'descend');      % target distribution (not necessary to be 'descend'
% rho0 = ones(1,N);
% rho0 = rho0/sum(rho0);
rho0 = rand(1,N);
rho0 = rho0/sum(rho0);

tspan = [0,100];
% deltaT = 0.025;
deltaT = 1e-3;
TotIt = tspan(2)/deltaT;
halftime = tspan(2)/2;
t = [tspan(1):deltaT:tspan(2)]';
samplesize = 1e4;

mode = 'Print';


%% Create a Q-matrix and run MH-ode iteration 

Q = QMH(pai);
edge = edgeQ(pai);
Q = Q.*edge;
Qrow = RowSumZero(Q)

rhoMHode = Iter_MHjump(pai, rho0, Qrow, tspan, deltaT, samplesize, mode);




