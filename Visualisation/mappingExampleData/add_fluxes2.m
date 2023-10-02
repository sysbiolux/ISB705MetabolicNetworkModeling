clearvars -except solverOK, clc, close all
load consistent_model.mat
load('allstatsbySubsystems.mat')
load('annotation_cyto.mat')

annotationscytoscape(1,:) = [];

dico=[consistent_model.rxnNames,consistent_model.rxns, strcat('R_',consistent_model.rxns), consistent_model.subSystems]
dico(:,3)=strrep(dico(:,3),'[','__91__')
dico(:,3)=strrep(dico(:,3),']','__93__')
dico(ismember(dico(:,2),'PCHOLSTE_HSABCt'),3)= cellstr('R_PCHOLSTE_HSABCt')
dico(ismember(dico(:,2),'PCHOLSTE_HSt1e'),3)= cellstr('R_PCHOLSTE_HSt1e')
dico(ismember(dico(:,2),'PCHOLPALME-HSABCt'),3)= cellstr('R_PCHOLPALME__45__HSABCt')
dico(ismember(dico(:,2),'PCHOLPALME-HSt1e'),3)= cellstr('R_PCHOLPALME__45__HSt1e')

dico(ismember(dico(:,2),'ALA-DTDe'),3)= cellstr('R_ALA__45__DTDe')
missing=setdiff(annotationscytoscape.sbmlid,dico(:,3))
save('dico4cytoscape','dico')

dico=array2table(dico);
[II, IIA, IIB]=intersect(allstatsbySubsystems.rxns,strrep(setdiff(dico.dico2,allstatsbySubsystems.rxns),'_r',''))
tmp=strrep(setdiff(dico.dico2,allstatsbySubsystems.rxns),'_r','');
allstatsbySubsystems.rxns(IIA) = strcat(tmp,'_r');
[I,IA, IB]=intersect(allstatsbySubsystems.rxns,dico.dico2);
dico(IB,end+1 : end + size(allstatsbySubsystems,2) )=allstatsbySubsystems(IA,:);
dico.Properties.VariableNames(9:end)=allstatsbySubsystems.Properties.VariableNames(5:end);
writetable(dico(:,[1,9]))
