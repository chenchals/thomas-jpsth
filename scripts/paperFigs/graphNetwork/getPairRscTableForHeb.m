function [hebTable] = getPairRscTableForHeb()
%CREATEPAIRRSCTABLEFORHEB 
% Uses cols below for significance flags
%     rhoCol = 'rhoRaw_150ms';
%     pvalCol = 'pvalRaw_150ms';
%% Load static RSC matrix for all pairs of SEF units and create file for Hierarchical Edge Bundling plot in R
   % Use columns:
    rhoCol = 'rhoRaw_150ms';
    pvalCol = 'pvalRaw_150ms';
    pvalThreshold = 0.01;

    fnUnitInfo = 'dataProcessed/dataset/dataNeurophys_SAT.mat';% for unit info
    fnRho = 'dataProcessed/analysis/spkCorr/summary/spkCorrAllPairsStaticRhoPval.mat';
        % Load data
    temp = load(fnRho,'-regexp'); % load cross areas and self areas
    temp.spkCorrColumnDefs=[];
    allSpkCorr = table();
    fns = fieldnames(temp);
    for jj = 1:numel(fns)
        allSpkCorr = [allSpkCorr;temp.(fns{jj})]; %#ok<*AGROW>
    end
    unitInfo = load(fnUnitInfo,'unitInfo');
    unitInfo = unitInfo.unitInfo;
  
    % add functional types for X and Y units: unitInfo colNames are
    funcTypeCols = {'visGrade','visType','moveGrade','errGrade','rewGrade','poorIso'};
    funcTypeColsRegex = strcat('^',funcTypeCols,'$');
    allSpkCorr = innerjoin(allSpkCorr,unitInfo,'LeftKeys','X_unitNum',...
        'RightKeys','unitNum','RightVariables',funcTypeCols);
    allSpkCorr.Properties.VariableNames = regexprep(allSpkCorr.Properties.VariableNames,funcTypeColsRegex,strcat('X_',funcTypeCols));
    allSpkCorr = innerjoin(allSpkCorr,unitInfo,'LeftKeys','Y_unitNum',...
        'RightKeys','unitNum','RightVariables',funcTypeCols);
    allSpkCorr.Properties.VariableNames = regexprep(allSpkCorr.Properties.VariableNames,funcTypeColsRegex,strcat('Y_',funcTypeCols));
    
  
    % add monkey, sessNum, and sess
    allSpkCorr = innerjoin(allSpkCorr,unitInfo,'LeftKeys','X_unitNum',...
        'RightKeys','unitNum','RightVariables',{'monkey','sessNum','sess'});
    
    %clearvars temp fn* jj unitInfo
    allSpkCorr.outcome = regexprep(allSpkCorr.condition,'Fast|Accurate','');
    allSpkCorr.satCondition = regexprep(allSpkCorr.condition,'Correct|Error.*','');
    allSpkCorr.epoch = allSpkCorr.alignedName;
    allSpkCorr.sameAreaPair(strcmp(allSpkCorr.X_area,allSpkCorr.Y_area)) = 1;
    allSpkCorr.sameChannelPair(allSpkCorr.XY_Dist == 0) = 1;
    allSpkCorr.pairCount = ones(size(allSpkCorr,1),1);

    hebTable = allSpkCorr;
    signifPval = hebTable.(pvalCol) <= pvalThreshold;
    plusRho = hebTable.(rhoCol) >= 0;
    minusRho = hebTable.(rhoCol) < 0;
    
    hebTable.signifPlusRho(plusRho & signifPval) = 1 ;
    hebTable.signifMinusRho(minusRho & signifPval) = 1 ;
    hebTable.nonSignifRho(~signifPval) = 1;
    
    % for Hierarchical Edge Bundling plot in R
    % writetable(tTbl,'rcode/spkCorrVals.csv')
end

