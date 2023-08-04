%driver script
clear all
addpath(genpath('C:\Program Files\IBM\ILOG\\CPLEX_Studio1210\cplex\matlab\x64_win64'))
changeCobraSolver('ibm_cplex')
model=readCbModel('C:\Users\maria.pacheco\Documents\GitHub\basicAnalysis\Model\Blautia_hydrogenotrophica_DSM_10507.xml');
[TableNumerics, Required4biomass, Results, minimalMedia, model]=InitialQualityControl(model, 'EX_biomass(e)');