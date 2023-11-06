clearvars -except solverOK, clc, close all
load consistent_model.mat
model_orig=consistent_model;
model=model_orig;

epsilon=1e-4

FBAsolution=optimizeCbModel(model,'max','zero')

%% inputs
% model=model_cons
% v=sol.v %table2array(res(:,3)); %data to visualize
model=model_orig
v=FBAsolution.v %table2array(res(:,3)); %data to visualize

cutoffFluxSum=eps %keep metabolites above this cutoff
cutoffRxn=eps %keep rxn above this cutoff

keepSubSystemsOnly=[];
% keepSubSystemsOnly={'Glycolysis/gluconeogenesis'}
% keepSubSystemsOnly={'Glycolysis/gluconeogenesis','Pentose phosphate pathway','Citric acid cycle'}

removeCofactorsFromFile=[];
% removeCofactorsFromFile='metsCofactors.txt'

keepAllRxnsFromMets=0 %add all rxns for included metabolites, even if not in subsystem

maxNodeSize=150 %scale nodeSize according to this value
log2NodeSize=1 %1: log2(v+2), 0: no log2 scaling

outputFile='metsFluxSum_log2'

%% calculate flux sum per metabolite
temp=repmat(v',size(model.S,1),1);
fluxes=model.S.*temp;
fluxSumP=full(sum((fluxes>0).*fluxes,2));
fluxSumN=full(sum((fluxes<0).*fluxes,2));
temp=[fluxSumP, fluxSumN];

% number of overall and active reactions per metabolite
groupCount=full(sum(model.S~=0,2));
groupCountFlux=full(sum(fluxes~=0,2));

% rename metabolites for IDARE: coa[c] to M_coa__91__c__93__
temp=strcat('M_', model.mets);
temp=strrep(temp,'[','__91__');
temp=strrep(temp,']','__93__');

% nodeSize and shape and organize in table
if log2NodeSize
    nodeSize=log2(fluxSumP+2);
else
    nodeSize=fluxSumP;
end
keepMaxValue=max(nodeSize);
nodeSize=nodeSize*maxNodeSize/keepMaxValue; %scale to maxNodeSize
G=table(model.mets,temp,groupCount,groupCountFlux,fluxSumP,fluxSumN,nodeSize);
G.shape=repmat(cellstr('Elipse'),size(G,1),1);
G(1:10,:)

% top hits (largest flusSum)
Gup = sortrows(G,4,'descend');
Gdn = sortrows(G,4,'ascend');
Gup(1:30,:)
% Gdn(1:10,:)

% add keep/remove flag (keep above fluxSum cutoff)
keepremove={};
for counter=1:numel(fluxSumP)
    if fluxSumP(counter)>=cutoffFluxSum
        keepremove=[keepremove; 'keep'];
    else
        keepremove=[keepremove; 'remove'];
    end
end
keepremove(1:10)
G.keepremove=keepremove;

writetable(G(:,[2,7]),'metsFluxSum_log2',"QuoteStrings",0,'WriteVariableNames',0)

%% reactions
name=strcat('R_', model.rxns);
% R_EX_pcholn203_hs__91__e__93__
name=strrep(name,'[','__91__');
name=strrep(name,']','__93__');
name=strrep(name,'-','__45__');

% % % dico(ismember(dico(:,2),'PCHOLSTE_HSABCt'),3)= cellstr('R_PCHOLSTE_HSABCt');
% % % dico(ismember(dico(:,2),'PCHOLSTE_HSt1e'),3)= cellstr('R_PCHOLSTE_HSt1e');

if log2NodeSize
    nodeSize=log2(abs(v)+2);
else
    nodeSize=abs(v);
end
nodeSize=nodeSize*maxNodeSize/keepMaxValue; %scale to maxNodeSize (of mets)
shape={};
for counter=1:numel(v)
    if v(counter)>=cutoffRxn
        shape=[shape; 'Triangle'];
    elseif v(counter)<=-cutoffRxn
        shape=[shape; 'V'];
    else
        shape=[shape; 'none'];
    end
end
shape(1:10)

keepremove={};
for counter=1:numel(v)
    if abs(v(counter))>=cutoffRxn
        keepremove=[keepremove; 'keep'];
    else
        keepremove=[keepremove; 'remove'];
    end
end
keepremove(1:10)

out1=G(:,[2,7,8,9]);
out2=table(name, nodeSize, shape, keepremove,'VariableNames',out1.Properties.VariableNames);
out1(1:10,:)
out2(1:10,:)
out=[out1; out2];

%% subSystems only / remove cofactors / keep all rxns from mets
toKeep=[];
if ~isempty(keepSubSystemsOnly)
    for counter=1:numel(keepSubSystemsOnly)
        subsystem=keepSubSystemsOnly(counter)
        temp=find(ismember(model.subSystems,subsystem));
        temp2=findMetsFromRxns(model,model.rxns(temp)); %mets
        % rxns above cutoff only
        vs=[];
        for counter2=1:numel(temp)
            %             temp2=find(ismember(model.rxns,temp(counter2)));
            vs=[vs; v(temp(counter2))];
        end
        %        vs
        temp(abs(vs)<cutoffRxn)=[];
        
        fs=[];
        for counter2=1:numel(temp2)
            temp3=find(ismember(model.mets,temp2(counter2)));
            fs=[fs; fluxSumP(temp3)];
        end
        % fs
        temp2(fs<cutoffFluxSum)=[];
        
        temp=strcat('R_',model.rxns(temp)); %rxns to keep
        temp=strrep(temp,'[','__91__');
        temp=strrep(temp,']','__93__');
        temp=strrep(temp,'-','__45__');
        
        temp2=strcat('M_', temp2); %mets to keep
        temp2=strrep(temp2,'[','__91__');
        temp2=strrep(temp2,']','__93__');
        
        toKeep=[toKeep; unique([temp;temp2])];
    end
    temp3=find(ismember(table2cell(out(:,1)),unique(toKeep)));
    out(temp3,4)=repmat({'keep'},numel(temp3),1);
    
    temp4=setdiff(1:size(out,1),temp3);
    out(temp4,4)=repmat({'remove'},numel(temp4),1);
end

if ~isempty(removeCofactorsFromFile) %remove cofactor metabolites
    fileID = fopen(removeCofactorsFromFile);
    C = textscan(fileID,'%s','Delimiter','\n')
    C=C{:}
    fclose(fileID);
    
    temp=find(ismember(model.metNames,C));
    temp=model.mets(temp);
    temp=strcat('M_', temp);
    temp=strrep(temp,'[','__91__');
    temp=strrep(temp,']','__93__');
    temp2=find(ismember(table2cell(out(:,1)),temp));
    out(temp2,4)=repmat({'remove'},numel(temp2),1);
end

if keepAllRxnsFromMets
    temp=find(ismember(table2cell(out(:,4)),'keep'));
    temp2=find(contains(table2cell(out(:,1)),'M_'));
    temp2(temp2>numel(model.mets))=[]; %no rxns here
    temp3=intersect(temp,temp2);
    temp4=cellstr(out{temp3,1});
    temp5=eraseBetween(temp4,1,2);
    temp5=strrep(temp5,'__91__','[');
    temp5=strrep(temp5,'__93__',']') %all mets
    
    temp=findRxnsFromMets(model,temp5);  %rxns to keep
    %above cutoff only
    vs=[];
    for counter=1:numel(temp)
        temp2=find(ismember(model.rxns,temp(counter)));
        vs=[vs; v(temp2)];
    end
    temp(abs(vs)<cutoffRxn)=[];
    
    temp=strcat('R_',temp);
    temp=strrep(temp,'[','__91__');
    temp=strrep(temp,']','__93__');
    temp=strrep(temp,'-','__45__');
    
    %     temp2=strcat('M_', temp2); %mets to keep
    %     temp2=strrep(temp2,'[','__91__');
    %     temp2=strrep(temp2,']','__93__');
    
    temp3=find(ismember(table2cell(out(:,1)),temp));
    out(temp3,4)=repmat({'keep'},numel(temp3),1);
end

% final stats
temp=find(ismember(table2cell(out(:,4)),'keep'));
temp2=find(contains(table2cell(out(:,1)),'M_'));
temp2(temp2>numel(model.mets))=[]; %no rxns here
temp3=intersect(temp,temp2);
disp('')
disp('Kept metabolites:')
numel(temp3)
temp2=find(contains(table2cell(out(:,1)),'R_'));
temp3=intersect(temp,temp2);
disp('')
disp('Kept reactions:')
numel(temp3)

writetable(out,outputFile,"QuoteStrings",0,'WriteVariableNames',0)
