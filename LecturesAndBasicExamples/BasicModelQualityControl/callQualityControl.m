%% driver script
clearvars -except solverOK, clc, close all
% initialize cplex if needed: 
% addpath(genpath('C:\Program Files\IBM\ILOG\\CPLEX_Studio1210\cplex\matlab\x64_win64'))
% changeCobraSolver('ibm_cplex')

%% read model 1 (adapt file path if necessary)
model=readCbModel('..\..\Model\Blautia_hydrogenotrophica_DSM_10507.xml')
% fastbox toolbox needed here (add to path):
[TableNumerics, Required4biomass, Results, minimalMedia, model, consistentModel]=InitialQualityControl(model, 'EX_biomass(e)');

%% read model 2 (adapt file path if necessary)
load('..\..\Model\brain_model.mat')
% fastbox toolbox needed here (add to path):
[TableNumerics2, Required4biomass2, Results2, minimalMedia2, model2, consistentModel2]=InitialQualityControl(brain_model, 'biomass_maintenance');
