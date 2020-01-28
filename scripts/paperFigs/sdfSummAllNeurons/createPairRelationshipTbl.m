% Create a tables with the following fields
% table1: cross area pairs: SEF-FEF, SEF-SC, and FEF-SC
% table2: same area pairs: SEF-SEF, FEF-FEF, and SC-SC  
% sessId|unitNum|unitArea|SatCondition|epoch|signifPlusRhoCount|...
%        signifMinusRhoCount|nonSignifCount|sameChannel
% Data that needs to be loaded are for all pairs from
% summary/spkCorrAllPairsStaticRhoPval.mat 
%% Load static RSC matrix for all pairs of SEF units
    fnUnitInfo = 'dataProcessed/dataset/dataNeurophys_SAT.mat';% for unit info
    fnRho = 'dataProcessed/analysis/spkCorr/summary/spkCorrAllPairsStaticRhoPval.mat';
        % Load data
    temp = load(fnRho,'-regexp'); % load cross areas and self areas
    allSpkCorr = table();
    fns = fieldnames(temp);
    for jj = 1:numel(fns)
        allSpkCorr = [allSpkCorr;temp.(fns{jj})]; %#ok<*AGROW>
    end
    unitInfo = load(fnUnitInfo,'unitInfo');
    unitInfo = unitInfo.unitInfo;
    % add monkey, sessNum, and sess
    allSpkCorr = innerjoin(allSpkCorr,unitInfo,'LeftKeys','X_unitNum',...
        'RightKeys','unitNum','RightVariables',{'monkey','sessNum','sess'});
    clearvars temp fn* jj unitInfo
    allSpkCorr.outcome = regexprep(allSpkCorr.condition,'Fast|Accurate','');
    allSpkCorr.satCondition = regexprep(allSpkCorr.condition,'Correct|Error.*','');
    allSpkCorr.epoch = allSpkCorr.alignedName;
    allSpkCorr.sameAreaPair(strcmp(allSpkCorr.X_area,allSpkCorr.Y_area)) = 1;
    allSpkCorr.sameChannelPair(allSpkCorr.XY_Dist == 0) = 1;
    allSpkCorr.pairCount = ones(size(allSpkCorr,1),1);
    %%
    % Use columns:
    rhoCol = 'rhoRaw_150ms';
    pvalCol = 'pvalRaw_150ms';
    pvalThreshold = 0.01;
    tTbl = allSpkCorr;
    signifPval = tTbl.(pvalCol) <= pvalThreshold;
    plusRho = tTbl.(rhoCol) >= 0;
    minusRho = tTbl.(rhoCol) < 0;
    
    
    tTbl.signifPlusRho(plusRho & signifPval) = 1 ;
    tTbl.signifMinusRho(minusRho & signifPval) = 1 ;
    tTbl.nonSignifRho(~signifPval) = 1;

    % seperate X_unit and Y_unit anddd aggregate as unitNum column table
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
    sameAreaPairs = 0;
    temp = useTbl(useTbl.sameAreaPair==sameAreaPairs,:);
    crossAreaTbl = grpstats(temp, ...
        {'monkey','sess','sessNum','unitNum','unitArea','satCondition','outcome','epoch'},{'sum'});
    
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
    
    
    
    
 
    
    
    
    
    

