clear; close all force; clc;
addpath('../')
load([pwd '\evap53pts.mat'])
ps.doeBuild = doeBuild;
% dace parameters
ps.regModel = @regpoly1;
ps.corModel = @corrgauss;

% Define domain bounds
ps.lb = [8.5 0  102 0];    % for evaporator
ps.ub = [20 100 400 400]; % for evaporator

% Define column index of design, objective function and constraint values
ps.designIndex = 4:7; % designIndex = 4:7 for evaporator
ps.fobjIndex = 20;     % 20 for evaporator
ps.constIndex = [11 14]; % 11 14 for evaporator

% define constraint limit and number of constraints
ps.constLimit = [35.5 40 80]; % [35.5 40 80] for evaporator

% Specify constraint tolerance. g(x) <= conTol
ps.conTol = 1e-6;

% define constraint type
ps.constType = [1 2];

ps.tol1 = 1e-4;
ps.tol2 = 1e-6;

ps.tolContraction = 1e-4;

ps.firstContractionFactor = 0.8;
ps.secondContractionFactor = 0.4;

ps.MaxFunEvals = 100*length(ps.lb);

ps.penaltyFactor = 1e3;

ps.nlp_solver = 'ipopt';

caballero(@evaporatorDOE, ps)