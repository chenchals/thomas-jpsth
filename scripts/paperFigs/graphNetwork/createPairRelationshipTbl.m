% Create a tables with the following fields
% table1: cross area pairs: SEF-FEF, SEF-SC, and FEF-SC
% table2: same area pairs: SEF-SEF, FEF-FEF, and SC-SC  
% sessId|unitNum|unitArea|SatCondition|epoch|signifPlusRhoCount|...
%        signifMinusRhoCount|nonSignifCount|sameChannel
% Data that needs to be loaded are for all pairs from
% summary/spkCorrAllPairsStaticRhoPval.mat 
% create file for Hierarchical Edge Bundling plot in R 

%% Load static RSC matrix for all pairs of SEF units and create file for Hierarchical Edge Bundling plot in R

    tTbl = getPairRscTableForHeb();    
    % for Hierarchical Edge Bundling plot in R
    % writetable(tTbl,'rcode/spkCorrVals.csv')
    
%% seperate X_unit and Y_unit anddd aggregate as unitNum column table
    % table size will become 2x
    tTbl1 = tTbl;
    tTbl1.unitNum = tTbl1.X_unitNum;
    tTbl1.unitArea = tTbl.X_area;
    tTbl2 = tTbl;
    tTbl2.unitNum = tTbl2.Y_unitNum;
    tTbl2.unitArea = tTbl2.Y_area;
    useTbl = [tTbl1;tTbl2];
        useCols = {
        'monkey'
        'sessNum'
        'sess'
        'unitNum'
        'unitArea'
        'satCondition'
        'outcome'
        'epoch'
        'sameAreaPair'
        'sameChannelPair'
        'signifPlusRho'
        'signifMinusRho'
        'nonSignifRho'
        'pairCount'
        };
    useTbl = useTbl(:,useCols);
    %%
    % group cross area pairs
    fx_renmedCols= @(tbl) regexprep(tbl.Properties.VariableNames, 'sum_','');
    sameAreaPairs = 0;
    temp = useTbl(useTbl.sameAreaPair==sameAreaPairs,:);
    crossAreaTbl = grpstats(temp, ...
        {'monkey','sess','sessNum','unitNum','unitArea','satCondition','outcome','epoch'},{'sum'});
    %crossAreaTbl.Properties.VariableNames = regexprep(crossAreaTbl.Properties.VariableNames, 'sum_', '');
    crossAreaTbl.Properties.VariableNames = fx_renmedCols(crossAreaTbl);
    %% Plot the counts of neurons with (+) significant, (-) significant, and non-significant correlation,
    %% by area (SEF, FEF, SC), for the PostSaccade epoch during Correct trials.
    
    tmpTbl = crossAreaTbl;
%     tmpTbl.Properties.RowNames = {}; %clear out row names in preparation for groupStats
    
    epoch = 'PostSaccade';
    outcome = 'Correct';
    tmpTbl = tmpTbl(ismember(tmpTbl.epoch, epoch) & ismember(tmpTbl.outcome, outcome), :);
    
    plotCols = {
      'unitArea'
      'satCondition'
      'signifPlusRho'
      'signifMinusRho'
      'nonSignifRho'
      'pairCount'};
    %t = tmpTbl(:,plotCols);
    
    outTbl = grpstats(tmpTbl(:,plotCols), {'satCondition','unitArea'}, {'sum'});
    outTbl.Properties.VariableNames = fx_renmedCols(outTbl);
    
    %% Function to spot check...
    function [] = checkNotYet(inFinalCountsTbl,unitNum,outcome,epoch)
        % check
    Z=inFinalCountsTbl(ismember(inFinalCountsTbl.epoch,epoch) & ismember(inFinalCountsTbl.outcome,outcome),:);
    Z=sortrows(Z,'unitArea','descend');
    Z = Z(Z.unitNum==unitNum,:);
    % check
    Z2 = temp(temp.unitNum==5 & ismember(temp.outcome,'Correct') & ismember(temp.epoch,'PostSaccade'),:);
    Z2Fast = Z2(ismember(Z2.satCondition,'Fast'),:);
    Z2Accurate = Z2(ismember(Z2.satCondition,'Accurate'),:);
    
    
    Z3 = tTbl((tTbl.X_unitNum==5 | tTbl.Y_unitNum==5) ...
              & tTbl.sameAreaPair == 0 ...
              & ismember(tTbl.outcome,'Correct') ...
              & ismember(tTbl.epoch,'PostSaccade'),:);

    
    end
    
    
    
    
 
    
    
    
    
    

