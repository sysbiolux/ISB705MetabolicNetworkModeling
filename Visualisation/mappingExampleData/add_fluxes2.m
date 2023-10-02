clearvars -except solverOK, clc, close all
load consistent_model.mat
load('allstatsbySubsystems.mat')
load('annotation_cyto.mat')

annotationscytoscape(1,:) = [];

dico=[consistent_model.rxnNames,consistent_model.rxns, strcat('R_',consistent_model.rxns), consistent_model.subSystems];
dico(:,3)=strrep(dico(:,3),'[','__91__');
dico(:,3)=strrep(dico(:,3),']','__93__');
dico(ismember(dico(:,2),'PCHOLSTE_HSABCt'),3)= cellstr('R_PCHOLSTE_HSABCt');
dico(ismember(dico(:,2),'PCHOLSTE_HSt1e'),3)= cellstr('R_PCHOLSTE_HSt1e');
dico(ismember(dico(:,2),'PCHOLPALME-HSABCt'),3)= cellstr('R_PCHOLPALME__45__HSABCt');
dico(ismember(dico(:,2),'PCHOLPALME-HSt1e'),3)= cellstr('R_PCHOLPALME__45__HSt1e');

dico(ismember(dico(:,2),'ALA-DTDe'),3)= cellstr('R_ALA__45__DTDe');
missing=setdiff(annotationscytoscape.sbmlid,dico(:,3));
save('dico4cytoscape','dico')

dico=array2table(dico);
[II, IIA, IIB]=intersect(allstatsbySubsystems.rxns,strrep(setdiff(dico.dico2,allstatsbySubsystems.rxns),'_r',''))
tmp=strrep(setdiff(dico.dico2,allstatsbySubsystems.rxns),'_r','');
allstatsbySubsystems.rxns(IIA) = strcat(tmp,'_r');
[I,IA, IB]=intersect(allstatsbySubsystems.rxns,dico.dico2);
dico(IB,end+1 : end + size(allstatsbySubsystems,2) )=allstatsbySubsystems(IA,:);
dico.Properties.VariableNames(9:end)=allstatsbySubsystems.Properties.VariableNames(5:end);
writetable(dico(:,[1,9]))

%%
res=dico(:,[1,9]);
size(res,1)
temp=(abs(table2array(res(:,2)))<10);
sum(temp)
toRemove=dico{temp,1};
writecell(toRemove,'rxnsToRemove.txt',"QuoteStrings",0)

%% use similarity between metabolites (only excact chemical similarity)
model=consistent_model;
% metScores = compareAllMets_TS(consistent_model,consistent_model);
load metscores
temp=find(contains(model.metNames,'Glucose'))
model.mets(temp)
model.metNames(temp)

% glc_D[c]
temp=find(contains(model.mets,'glc_D[c]'))
data=metScores(temp,:);
temp2=(data>0.001);
sum(temp2)
model.metNames(temp2)

%% use similarity between metabolites (approx)
include={'C6','C5','C4','C3'};
match=zeros(size(model.mets));
for counter=1:numel(include)
    name=include(counter)
    temp = regexpi(model.metFormulas, name);
    temp = ~cellfun(@isempty, temp);
    sum(temp)
    match=match+temp;
end
match=(match>0);
sum(match)
% model.metNames(find(match))
[rxnList, rxnFormulaList] = findRxnsFromMets(model, model.mets(match));

% combine with flux cutoff
res=dico(:,[1,2,9]);
size(res,1)
temp=(abs(table2array(res(:,3)))<10);
sum(temp)
res(temp,:)=[];

[C,IA]=setdiff(res{:,2},rxnList);
res(IA,:)=[];

[C,IA]=setdiff(model.rxns, res{:,2});
toRemove2=dico{IA,1};
toRemove2=[toRemove2; model.metNames(find(match==0))]; %add non-matching metabolites
writecell(toRemove2,'ToRemove2.txt',"QuoteStrings",0)
