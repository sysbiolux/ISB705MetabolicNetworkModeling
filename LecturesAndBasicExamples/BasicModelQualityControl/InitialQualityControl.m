function[TableNumerics, Required4biomass, Results, minimalMedia, model, ConsistentModel]=InitialQualityControl(model, objective)
% (c) Maria Pires Pacheco and Thomas Sauter (University of Luxembourg)
[model,ConsistentModel, ~,~ , ~]=checkModel(model);
fields2check={'mets','rxns','genes'}
for i=1:numel(fields2check)
    if numel(unique(model.(fields2check{i}))) == numel(model.(fields2check{i}))
    else
        [~,idxu,idxc] = unique(model.(fields2check{i}));
        % count unique values (use histc in <=R2014b)
        [count, ~, idxcount] = histcounts(idxc,numel(idxu));
        % Where is greater than one occurence
        idxkeep = count(idxcount)>1;
        item=model.(fields2check{i});
        to_remove=item(idxkeep);
        if i==2
            model=removeRxns(model,to_remove(2:end));
        end
        
    end
    model=buildRxnGeneMat(model);
end
model.subSystems(cellfun('isempty',model.subSystems))=cellstr('');
if isfield(model,'subSystems')
    if iscellstr( model.subSystems)
    elseif iscell( model.subSystems)
        if size(model.subSystems,2)==1
            try tmp=vertcat(model.subSystems{:});
                
                if numel(tmp)==numel(model.subSystems)
                    model.subSystems=tmp;
                else
                    for i=1:numel(model.subSystems)
                        tmp=model.subSystems{i,1}{1};
                        if ischar(tmp)
                            tmp=cellstr(tmp);
                        end
                        model.subSystems(i)=tmp;
                    end
                end
            catch
                
                for i=1:numel(model.subSystems)
                    disp('Warning reactions are assigned to multiple pathways')
                    tmp=model.subSystems{i};
                    if ischar(tmp)
                        tmp=cellstr(tmp);
                    elseif iscell(tmp)
                        tmp=tmp(1);
                    end
                    model.subSystems(i)=tmp;
                end
            end
        end
    elseif ischar(model.subSystems)
    else
        disp('Warning there is no subSystems')
    end
    clear i count
    save
    model.subSystems(cellfun('isempty',model.subSystems))=cellstr('');
    TableNumerics=table(numel(model.mets),numel(model.rxns),numel(model.genes), numel(unique(model.subSystems)));
    TableNumerics.Properties.VariableNames={'mets', 'rxns', 'genes', 'pathways'};
    
else
    TableNumerics=table(numel(model.mets),numel(model.rxns),numel(model.genes));
    TableNumerics.Properties.VariableNames={'mets', 'rxns', 'genes'};
end
disp('Numerics:')
disp(TableNumerics)
closedRxns=model.rxns(model.lb==0 & model.ub==0);
disp('closed reactions');
disp(closedRxns)
Results.closed=closedRxns;

bounds_span=min(abs(model.ub-model.lb));
disp(strcat( 'bounds span: ',num2str(bounds_span)))
Results.bounds_span=bounds_span;


disp(strcat('number of objectives:', num2str(numel(find(model.c)))));
if sum(model.c ~=0 )>0
    model.rxns(model.c==1)
end

%% flux consistency

disp(strcat('consistency rate:', num2str(numel(ConsistentModel.rxns)/numel(model.rxns))));

blocked=setdiff(model.rxns,ConsistentModel.rxns);
disp(strcat('blocked reactions:',num2str(numel(blocked))));
Results.blocked=blocked;
%% growth rate
ConsistentModel = changeObjective(ConsistentModel,objective);
sol=optimizeCbModel(ConsistentModel,'max','zero'); %minimal solutio
sum(sol.x~=0) %number of non-zero fluxes
if ~isempty(sol.f)
if isfield(ConsistentModel, 'subSystems')
    temp=(model.subSystems(sol.v~=0));
    [uc, ~, idc] = unique( temp ) ;
    counts = accumarray( idc, ones(size(idc)) ) ;
    Required4biomass=table(uc,counts); %involved subSystems
else
    Required4biomass=[];
end
save
exR=find(findExcRxns(ConsistentModel))
[~,match]=find(ConsistentModel.S(:,exR)<0)
[~,match2]=find(ConsistentModel.S(:,exR)>0)

% minimal medium
minimalMedia=intersect(ConsistentModel.rxns(sol.v<0),ConsistentModel.rxns(exR(match)));
minimalMedia2=intersect(ConsistentModel.rxns(sol.v>0),ConsistentModel.rxns(exR(match2)));
minimalMedia=union(minimalMedia, minimalMedia2);
%secretion
secretion=intersect(ConsistentModel.rxns(sol.v>0),ConsistentModel.rxns(exR(match)));
secretion2=intersect(ConsistentModel.rxns(sol.v<0),ConsistentModel.rxns(exR(match2)));
secretion=union(secretion, secretion2);


Flux=table(ConsistentModel.rxns,sol.v);
Results.flux=Flux;
Results.secretion=secretion;
end
end