%% Building of context-specific models using rFASTCORMICS (Pacheco et al, 2019)
sampleModel=1; consensusModel=1;

if consensusModel==1  && sampleModel==1
    name ={'consensus', 'sample'};
elseif  consensusModel==1  && sampleModel==0
    name ={'consensus'};
elseif  consensusModel==0 && sampleModel==1
    name ={'sample'};
end
ModelMatrix=struct();
path=('C:\Users\maria.pacheco\Documents\');
%optional
set(0,'defaultTextInterpreter','none');
set(groot,'defaulttextinterpreter','none');
set(groot, 'DefaultAxesTickLabelInterpreter', 'none');
set(groot, 'defaultLegendInterpreter','none');
%Uncomment if an error related to eval occurs
%feature astheightlimit 2000
% color code used for the heat maps
altcolor= [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
    255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0]/255; %shorter 10% = 1 bar

delete clone*.log %delete some files generated by cplex
%% 3.1 Expression data input
% Import the gene expression data into Matlab.

data_cancer = readtable('fpkm_BRCA_cancer.txt', "ReadRowNames",true);
data_control = readtable('fpkm_BRCA_control.txt', "ReadRowNames",true);
data = [data_cancer, data_control];
fpkm = table2array(data);
rownames = data.Properties.RowNames;
colnames = data.Properties.VariableNames;
conditions=colnames;
conditions(contains(conditions,'BRCA_H'))=cellstr('Healthy');
conditions(contains(conditions,'BRCA_C'))=cellstr('Cancer');
uniqueconditions=unique(conditions);


%% 3.2. Context-specific model reconstruction
% human reconstruction Recon 3

load('Recon3Model.mat')
%% load medium - dictionnaries -

load mediumExample.mat
load dicorFASTCORMICS.mat

epsilon = 1e-4;
consensusProportion = 0.9;
%unnest the model.subSystems
if ischar(model.subSystems{1})
else
    model.subSystems=vertcat(model.subSystems{:});
end

unpenalizedSystems = {'Transport, endoplasmic reticular';
    'Transport, extracellular';
    'Transport, golgi apparatus';
    'Transport, mitochondrial';
    'Transport, peroxisomal';
    'Transport, lysosomal';
    'Transport, nuclear'};
unpenalized = model.rxns(ismember(model.subSystems,unpenalizedSystems));
optionalSettings.unpenalized = unpenalized;
optionalSettings.func = {'DM_atp_c_','biomass_reaction'}; % forced reactions
notMediumConstrained = 'EX_tag_hs(e)';% if no constrain is used please remove the field.
optionalSettings.notMediumConstrained = notMediumConstrained;
optionalSettings.medium = medium_example;% if no medium is used please remove the field.

biomassReactionName = {'biomass_reaction'};
printFigures=1;
changeCobraSolver('ibm_cplex');

%% Reconstruct the context-specific models
if sum(ismember(name,'consensus'))
    MatrixConsensus =zeros(numel(model.rxns), numel(uniqueconditions));

    for i=1:numel(uniqueconditions)
        match=ismember(conditions,uniqueconditions(i));
        [index]=rFastcormics4cobra(model, fpkm(:,match), rownames, colnames(match),dico, consensusProportion, epsilon, optionalSettings, biomassReactionName, printFigures, path);
        MatrixConsensus(index,i)=1;
        MatrixConsensus(index,i)=1;
    end
    ModelMatrix.consensus=MatrixConsensus;
    ModelMatrix.conditions=uniqueconditions;

end
ModelMatrix.consistentModel=model;

%% sample-specific models
if sum(ismember(name,'sample'))
    sampleModelsMatrix =zeros(numel(model.rxns), numel(colnames));
    for i = 1:numel(colnames) %for each sample
        [ sample_models] = rFastcormics4cobra(model, fpkm(:,i), rownames, colnames(i),dico, consensusProportion, epsilon, optionalSettings, biomassReactionName, printFigures);
        sampleModelsMatrix(sample_models,i)=1;
    end
    ModelMatrix.sample=sampleModelsMatrix;
    ModelMatrix.colnames=colnames;
end
save Reconstruction