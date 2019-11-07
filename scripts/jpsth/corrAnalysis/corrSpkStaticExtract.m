oFilename = 'dataProcessed/analysis/spkCorr/spkCorrAllPairsStatic.mat';
areaPairs = {
    'SEF-SEF' 
    'SEF-FEF'    
    'SEF-SC'     
    'FEF-FEF'    
    'FEF-SC'     
    'SC-SC'      
    };
corrMatDirs = strcat('dataProcessed/analysis/spkCorr/spkCorr_',areaPairs,'/mat');
corrDatFields = {
    'Pair_UID'        
    'X_monkey'        
    'X_sessNum'       
    'X_sess'          
    'X_unitNum'       
    'Y_unitNum'       
    'X_unit'          
    'Y_unit'          
    'X_area'          
    'Y_area'          
    'X_visGrade'      
    'Y_visGrade'      
    'X_visField'      
    'Y_visField'      
    'X_visType'       
    'Y_visType'       
    'X_moveGrade'     
    'Y_moveGrade'     
    'X_moveField'     
    'Y_moveField'     
    'X_errGrade'      
    'Y_errGrade' 
    'X_isErrGrade'      
    'Y_isErrGrade' 
    'X_errField'      
    'Y_errField'      
    'X_rewGrade'      
    'Y_rewGrade'      
    'X_isRewGrade'      
    'Y_isRewGrade'      
    'X_taskType'      
    'X_Hemi'          
    'Y_Hemi'          
    'X_Grid'          
    'Y_Grid'          
    'X_GridAP_ML'     
    'Y_GridAP_ML'     
    'X_Depth'         
    'Y_Depth'         
    'X_Depth0'        
    'Y_Depth0'        
    'X_newDepth'      
    'Y_newDepth'      
    'XY_Dist'         
    'isOnSameChannel' 
    % corr fields
    'condition'          
    'alignedName'        
    'alignedEvent'       
    'alignedTimeWin'     
    %'trialNosByCondition'
    'critRho10'          
    'critRho05'          
    'critRho01'          
    'critRho10_Z'        
    'critRho05_Z'        
    'critRho01_Z'        
    'rho_pval_static'    
    'rho_pval_win_Z'     
    'rho_pval_static_Z_baseline'  
    'rho_pval_static_Z_trial'  
    };
spkCorrStatic = struct();
tic
for d = 1:numel(corrMatDirs)
    areaPair = areaPairs{d};
    fprintf('Doing pairs for %s...',areaPair);
    areaPairField = strrep(areaPair,'-','_');
    srcFiles = dir([corrMatDirs{d},'/spkCorr_PAIR_*.mat']);
    srcFiles = strcat({srcFiles.folder}','/',{srcFiles.name}');
    nPairs = numel(srcFiles);
    outPairs = struct();
    parfor p = 1:nPairs
        pDat = table();
        srcFile = srcFiles{p};
        pDatStruct = load(srcFile);
        [~,fn,ext] = fileparts(srcFile);
        srcFile = [fn ext]; 
        nRows = size(pDatStruct.spikeCorr,1);
        temp = [repmat(pDatStruct.cellPairInfo,nRows,1) pDatStruct.spkCorr];
        tempFns = temp.Properties.VariableNames';
        pDat.srcFile = repmat({srcFile},nRows,1);
        pDat.pairAreas = repmat({areaPair},nRows,1);
        pDat.nTrials = cellfun(@(x) numel(x), temp.trialNosByCondition);
        temp.trialNosByCondition = [];
        fieldIdx = cellfun(@(x) find(strcmp(tempFns,x)),corrDatFields)
        pDat = [pDat,temp(:,corrDatFields)];
        pDat.XY_Dist = cell2mat(pDat.XY_Dist);
        % split rho, pval
        [pDat.rhoRaw,pDat.pvalRaw,pDat.signifRaw_05,pDat.signifRaw_01] = ...
            cellfun(@(x) deal(x(1),x(2),x(2)<=0.05,x(2)<=0.01),pDat.rho_pval_static);
        [pDat.rhoZBaseline,pDat.pvalZBaseline,pDat.signifZBaseline_05,pDat.signifZBaseline_01] = ...
            cellfun(@(x) deal(x(1),x(2),x(2)<=0.05,x(2)<=0.01),pDat.rho_pval_static_Z_baseline);
        [pDat.rhoZTrial,pDat.pvalZTrial,pDat.signifZTrial_05,pDat.signifZTrial_01] = ...
            cellfun(@(x) deal(x(1),x(2),x(2)<=0.05,x(2)<=0.01),pDat.rho_pval_static_Z_trial);
        outPairs(p).pDat = pDat;
    end
    spkCorrStatic.(areaPairField) = vertcat(outPairs.pDat); 
    fprintf('Done %.3f sec.\n',toc)
end
size(spkCorrStatic);
save(oFilename,'-v7.3','-struct','spkCorrStatic');
toc




