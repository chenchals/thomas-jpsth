oFilename = 'dataProcessed/analysis/spkCorr/rscSignalNoiseStatic.mat';
areaPairs = {
    'SEF-SEF' 
    'SEF-FEF'    
    'SEF-SC'     
    'FEF-FEF'    
    'FEF-SC'     
    'SC-SC'      
    'SEF-NSEFN'  
    'FEF-NSEFN'  
    'SC-NSEFN'   
    'NSEFN-NSEFN'
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
    'X_errField'      
    'Y_errField'      
    'X_rewGrade'      
    'Y_rewGrade'      
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
    'trialNosByCondition'
    'critRho10'          
    'critRho05'          
    'critRho01'          
    'critRho10_Z'        
    'critRho05_Z'        
    'critRho01_Z'        
    'rho_pval_static'    
    'rho_pval_win_Z'     
    'rho_pval_static_Z'  
    };
rscSignalNoise = struct();
for d = 1:numel(corrMatDirs)
    areaPair = areaPairs{d};
    srcFiles = dir([corrMatDirs{d},'/rscCorr_PAIR_*.mat']);
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
        temp = [repmat(pDatStruct.cellPairInfo,nRows,1) pDatStruct.spikeCorr];
        tempFns = temp.Properties.VariableNames';
        pDat.srcFile = repmat({srcFile},nRows,1);
        pDat.pairAreas = repmat({areaPair},nRows,1);
        pDat.nTrials = cellfun(@(x) numel(x), temp.trialNosByCondition);
        temp.trialNosByCondition = [];
        pDat = [pDat,...
                temp(:,ismember(tempFns,corrDatFields))];                
        outPairs(p).pDat= pDat;
    end
    rscSignalNoise(d).pDat = vertcat(outPairs.pDat);    
end
rscSignalNoise = vertcat(rscSignalNoise.pDat);
size(rscSignalNoise);
save(oFilename,'rscSignalNoise');




