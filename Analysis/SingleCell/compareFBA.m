function[ComparisonFlux2Modelskeep,res4tKeep, FluxPerClusterKeep,FluxPerPathwayKeep, FluxPerPathwayPerClusterKeep]=compareFBA(FBAsolutionKeep, Results_keep)

k=0;
%%
ComparisonFlux2Modelskeep=struct();
res4tKeep=struct();
FluxPerClusterKeep=struct();
FluxPerPathwayKeep=struct();
FluxPerPathwayPerClusterKeep=struct();

for counterA=1:numel(FBAsolutionKeep)
    for counterB=1:numel(FBAsolutionKeep)

        multicell_model=Results_keep(counterA).multi_cell_population;

        if counterA<counterB
          k=k+1;
            data1=FBAsolutionKeep(counterA).Analysis.Var2;
            model1=multicell_model;
            multicell_model=Results_keep(counterB).multi_cell_population;
            data2=FBAsolutionKeep(counterB).Analysis.Var2;
            model2=multicell_model;

            if iscell(model1.subSystems{1})
                model1.subSystems= vertcat(model1.subSystems{:});
                if numel(model1.subSystems)< numel(model1.rxns)
                    model1.subSystems(end +1: numel(model1.rxns))=cellstr('');
                end

            end
            subSystems1=model1.subSystems;
            subSystems1(cellfun('isempty',model1.subSystems))=cellstr('NaN');
            if iscell(model2.subSystems{1})

                model2.subSystems= vertcat(model2.subSystems{:});
                if numel(model2.subSystems)< numel(model2.rxns)
                    model2.subSystems(end +1: numel(model2.rxns))=cellstr('');
                end
            end
            subSystems2=model2.subSystems;
            subSystems2(cellfun('isempty',model2.subSystems))=cellstr('NaN');
            disp('loading & pre-processing done ...')

            %% consolidate into one table
            [allRxns,~,~] = union(model1.rxns,model2.rxns);
            allSubSystems=cell(size(allRxns));
            v1=zeros(size(allRxns));
            v2=zeros(size(allRxns));
            v1nan=nan(size(allRxns));
            v2nan=nan(size(allRxns));
            [~,IC,ID]=intersect(allRxns,model1.rxns);
            allSubSystems(IC)=subSystems1(ID);
            v1(IC)=data1(ID);
            v1nan(IC)=data1(ID);
            [~,IC,ID]=intersect(allRxns,model2.rxns);
            allSubSystems(IC)=subSystems2(ID);
            v2(IC)=data2(ID);
            v2nan(IC)=data2(ID);
            diff_v2_v1=v2-v1;
            abs_diff=abs(diff_v2_v1);

            LastTwo = cellfun(@(S) (S(end)), allRxns, 'Uniform', 0);
            LastTwo((ismember(LastTwo,']')))={'0'};
            ClusterNr=str2num(cell2mat(LastTwo));

            ComparisonFlux2Models=table(allRxns,ClusterNr,allSubSystems,v1nan,v2nan,v1,v2,diff_v2_v1,abs_diff);
            ComparisonFlux2Modelskeep(k).ComparisonFlux2Modelskeep=ComparisonFlux2Modelskeep;



            %% stats
            temp=allSubSystems;
            [us, ~, ids] = unique( temp ) ;
            counts = accumarray( ids, ones(size(ids)) ) ;
            FluxPerPathway=table(us,counts);
           FluxPerPathwayKeep(k).FluxPerPathway=FluxPerPathway;
            temp=cellstr(LastTwo);
            [uc, ~, idc] = unique( temp ) ;
            counts = accumarray( idc, ones(size(idc)) ) ;
            FluxPerCluster=table(uc,counts);
           FluxPerClusterKeep(k).FluxPerCluster=FluxPerCluster;

            %% number of reactions per subSystem and cluster
            res2=nan(numel(us),numel(uc));
            for counter=1:numel(us) %subSystems
                for counter2=1:numel(uc) %clusters
                    idx=(ids==counter).*(idc==counter2);
                    res2(counter,counter2)=sum(idx);
                end
            end
            FluxPerPathwayPerCluster=array2table(res2,'RowNames',us,'VariableNames',uc');
            FluxPerPathwayPerClusterKeep(k).FluxPerPathwayPerCluster=FluxPerPathwayPerCluster;
            % absolute flux difference per subSystem and cluster
            res3=nan(numel(us),numel(uc));
            for counter=1:numel(us) %subSystems
                for counter2=2:numel(uc) %clusters
                    idx=(ids==counter).*(idc==counter2);
                    temp=table2array(ComparisonFlux2Models(find(idx),'abs_diff'));
                    res3(counter,counter2)=sum(temp);
                end
            end
            %res3t=array2table(res3,'RowNames',us,'VariableNames',uc');

            %temp=sort(max(res3,[],2),'descend');
            %show=find(max(res3,[],2)>=temp(10)); %show top10

            % relative flux difference per subSystem and cluster
            res4=res3./res2;
            res4t=array2table(res4,'RowNames',us,'VariableNames',uc');

            %temp=sort(max(res4,[],2),'descend');
            %show=find(max(res4,[],2)>=temp(10)); %show top10
            %resTop10=res4t(show,:);
            res4tKeep(k).rest=res4t;
            


        end
    end
end