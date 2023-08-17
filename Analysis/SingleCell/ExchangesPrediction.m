function[AnalysisKeep,ExchangesKeep]=ExchangesPrediction(Results_keep, model, T, CellNumberPerCluster)
AnalysisKeep=struct();
ExchangesKeep=struct();

%%
%PD
%STEP 1
%load the inputs in (model + model cluster names) and create a string with
%the cluster names + run FVA
for clusterrouned=1:numel(Results_keep)
modelNr=clusterrouned;


% multiCellModel=Results_keep(modelNr).multiCellModel
multiCellModel=Results_keep(modelNr).multi_cell_population;
LastTwo = cellfun(@(S) S(end-1:end), multiCellModel.rxns, 'Uniform', 0);
clusters=cell(size(LastTwo));
match=contains(LastTwo,'_');
clusters(match) = cellfun(@(S) S(end), LastTwo(match), 'Uniform', 0);

LastThree = cellfun(@(S) S(end-2:end), multiCellModel.rxns, 'Uniform', 0);
match=~contains(LastTwo,'_') & ~contains(LastTwo,'[e]');
clusters(match) = cellfun(@(S) S(end-1), LastThree(match), 'Uniform', 0);




[~,IA,IB]=intersect(multiCellModel.rxns,T.medium_rxns_keep);
multiCellModel.ub(IA)=T.met_Conc_mM(IB);
multiCellModel.ub(ismember(multiCellModel.rxns,'Ex_lnlc[e]'))=1;
multiCellModel.ub(ismember(multiCellModel.rxns,'Ex_o2s[e]'))=0; %superoxide
multiCellModel.ub(ismember(multiCellModel.rxns,'Ex_o2[e]'))=21.253*2; %o2, Glc x 2
multiCellModel.ub(ismember(multiCellModel.rxns,'Ex_co2[e]'))=21.253*2; %co2, Glc x 2
clear ans clear IA IB 

%setting a 10% lower bounds for each cluster's biomass reaction (10% of the
%max biomass number obtained when we fully optimise for the specific
%biomass reaction)
biomass=multiCellModel.rxns(contains(multiCellModel.rxns,'biomass'));
%for all metabolites that contains [u], R are the reactions that involve
%that metabolite -->

for i=1: size(CellNumberPerCluster,2)
    ordre = zeros(1,size(CellNumberPerCluster,2));
    if i==i
        ordre(i) = 1;
    end
    multiCellModel.c(contains(multiCellModel.rxns,'biomassmaintenance')) = ordre;
    %Flux Balance Analysis (FBA) is used to find the max optimised value if the
    %model was optimising only for this specific biomass reaction
    FBAsolMax1cluster=optimizeCbModel(multiCellModel);
    tenpercent = 0.1*FBAsolMax1cluster.f;
    %set lower bound
    multiCellModel.lb(contains(multiCellModel.rxns,biomass(i))) = tenpercent;
end
clear i ordre
TotalCells=sum(table2array(CellNumberPerCluster(clusterrouned,:)));
BiomassFix=table2array(CellNumberPerCluster(clusterrouned,:))./repmat(TotalCells,size(CellNumberPerCluster(clusterrouned,:),1),1);
multiCellModel.c(contains(multiCellModel.rxns,'biomassmaintenance')) = BiomassFix;
FBAweighted = optimizeCbModel(multiCellModel,'max','zero');
indx_non_zero_flux=FBAweighted.v~=0;
clear BiomassFix biomass FBAsolMax1cluster tenpercent TotalCells

ReactionsNonZerFluxFBA = multiCellModel.rxns(indx_non_zero_flux);
Flux_non_zero_FBA = FBAweighted.v(indx_non_zero_flux);
TableFlux=table(ReactionsNonZerFluxFBA, Flux_non_zero_FBA);
clear Flux_non_zero_FBA indx_non_zero_flux
%--> finding the U reactions and searching for the index so that we  can find
%the metabolites present in it with the S matrix
UMetsID=find(contains(multiCellModel.mets,'[u]'));


%--> list of metabolites  that take part in  U reactions and hav e a flux
%(those represent [u]metabolites and other metabolites present in U
%reactions)

%%

ClusterExchanges = zeros(length(UMetsID),size(CellNumberPerCluster,1)+1);
% ClusterExchanges.Properties.RowNames =  list_of_U_metabolites_that_have_flux(:,1);
% ClusterExchanges.Properties.VariableNames = list_U_reactions_FBAflux(:,1);

for i = 1:length(UMetsID)
    
    [~,Rindx]=find(multiCellModel.S(UMetsID(i),:));
    
    tmp=zeros(numel(Rindx),size(CellNumberPerCluster,1)+1);
    %to note down the flux of the reactions for that specific metabolite
    %decide if the metabolite is consumed or produced depending on which
    %side of the reaction it is (S matrix value) and depending on its flux
    % --> THIS  CODE HAS BEEN CORRECTED THROUGH THE TROUBLE SHOOTING PART
    % (SEE DOWN BELOW)
    for v=1:length(Rindx)
        %right and positive flux --> produced
        [~,~,IB]=intersect(multiCellModel.rxns(Rindx(v)),TableFlux.ReactionsNonZerFluxFBA);
        if ~isempty(IB)
            
            c=(clusters(Rindx(v)));
            if ~strcmp(c,'e')
                c=str2double(c{1});
            else
                c=size(CellNumberPerCluster,1)+1;
            end
            if multiCellModel.S(UMetsID(i),Rindx(v)) > 0 && TableFlux.Flux_non_zero_FBA(IB) > 0
                tmp(v,c)= table2array(TableFlux(ismember(TableFlux.ReactionsNonZerFluxFBA,multiCellModel.rxns(Rindx(v))),2));
                %left and positive flux --> consumed
            elseif multiCellModel.S(UMetsID(i),Rindx(v)) < 0 && TableFlux.Flux_non_zero_FBA(IB) > 0
                tmp(v,c)= -table2array(TableFlux(ismember(TableFlux.ReactionsNonZerFluxFBA,multiCellModel.rxns(Rindx(v))),2));
                %right and negative flux --> consumed
            elseif multiCellModel.S(UMetsID(i),Rindx(v)) > 0 && TableFlux.Flux_non_zero_FBA(IB) < 0
                tmp(v,c)= table2array(TableFlux(ismember(TableFlux.ReactionsNonZerFluxFBA,multiCellModel.rxns(Rindx(v))),2));
                %left and negative flux --> produced
            elseif multiCellModel.S(UMetsID(i),Rindx(v)) < 0 && TableFlux.Flux_non_zero_FBA(IB) < 0
                tmp(v,c)= -table2array(TableFlux(ismember(TableFlux.ReactionsNonZerFluxFBA,multiCellModel.rxns(Rindx(v))),2));
            end
            
            
            ClusterExchanges(i,:)=sum(tmp,1);
        end
    end
end
ClusterExchanges=array2table(ClusterExchanges);
ClusterExchanges.Properties.VariableNames=[CellNumberPerCluster.Row;cellstr('e')];
tmp=strrep(multiCellModel.mets(UMetsID),'[u]','[e]');
ClusterExchanges.Row=tmp;

[~,Ia,Ib]=(intersect(model.mets,tmp));
ClusterExchanges.Row(Ib)=model.metNames(Ia) ;

%cleaning the matrix from 0 fluxes or irrelevant low fluxes otherwise the
%plotting is too complicated and you can't see anything

Analysis=table(multiCellModel.rxns,FBAweighted.v, FBAweighted.x);
AnalysisKeep(clusterrouned).Analysis=Analysis;
ExchangesKeep(clusterrouned).Exchanges=ClusterExchanges;


end
end

