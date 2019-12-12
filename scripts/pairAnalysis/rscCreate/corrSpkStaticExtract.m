% Extract spike count correlation data for pairs of units tecorded form
% same session for same as well as different areas
% see also CREATECORRSPKBATCH

%%oFilename = 'dataProcessed/analysis/11-18-2019/spkCorr/summary/spkCorrAllPairsStaticNew.mat';
areaPairs = {
    'SEF-SEF' 
    'SEF-FEF'    
    'SEF-SC'     
    'FEF-FEF'    
    'FEF-SC'     
    'SC-SC'      
    };
corrMatDirs = strcat('dataProcessed/analysis/11-18-2019/spkCorr/spkCorr_',areaPairs,'/mat');
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
    'rho_pval_win_50ms'     
    'rho_pval_static_50ms'    
    'rho_pval_static_Z_baseline_50ms'  
    'rho_pval_static_Z_trial_50ms' 
    'rho_pval_win_150ms'     
    'rho_pval_static_150ms'    
    'rho_pval_static_Z_baseline_150ms'  
    'rho_pval_static_Z_trial_150ms' 
    'rho_pval_win_200ms'     
    'rho_pval_static_200ms'    
    'rho_pval_static_Z_baseline_200ms'  
    'rho_pval_static_Z_trial_200ms' 
    'xWaves'
    'yWaves'
    'xWaveWidths'
    'yWaveWidths'
    };

t = regexp(corrDatFields,'rho_pval_static_(\d*)ms$','tokens');
staticWinSizes = sort(cellfun(@(x) str2double(x{1}),t(~cellfun(@isempty,t))));

cols2Remove = {
    'xWaves'
    'yWaves'
    'xWaveWidths'
    'yWaveWidths'
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
        nRows = size(pDatStruct.spkCorr,1);
        temp = [repmat(pDatStruct.cellPairInfo,nRows,1) pDatStruct.spkCorr];
        tempFns = temp.Properties.VariableNames';
        pDat.srcFile = repmat({srcFile},nRows,1);
        pDat.pairAreas = repmat({areaPair},nRows,1);
        pDat.nTrials = cellfun(@(x) numel(x), temp.trialNosByCondition);
        temp.trialNosByCondition = [];
        fieldIdx = cell2mat(cellfun(@(x) find(strcmp(tempFns,x)),corrDatFields,'UniformOutput',false));
        pDat = [pDat,temp(:,corrDatFields)];
        pDat.XY_Dist = cell2mat(pDat.XY_Dist);
        %% split rho, pval
        for sW = 1:numel(staticWinSizes)
            staticWin = staticWinSizes(sW);
            swSuffix = num2str(staticWin,'_%dms');
            % Use rho/pval from raw counts            
            swRhoPval = pDat.(['rho_pval_static' swSuffix]);
            [rho,pval,sig05,sig01] = cellfun(@(x) deal(x(1),x(2),x(2)<=0.05,x(2)<=0.01),swRhoPval);            
            [pDat.(['rhoRaw' swSuffix]),...
            pDat.(['pvalRaw' swSuffix]),...
            pDat.(['signifRaw_05' swSuffix]),...
            pDat.(['signifRaw_01' swSuffix])] = deal(rho,pval,sig05,sig01);
            % Use rhp/pval from Z-scored (with mean and std from) Baseline period
            swRhoPval = pDat.(['rho_pval_static_Z_baseline' swSuffix]);
            [rho,pval,sig05,sig01] = cellfun(@(x) deal(x(1),x(2),x(2)<=0.05,x(2)<=0.01),swRhoPval);            
            [pDat.(['rhoZBaseline' swSuffix]),...
            pDat.(['pvalZBaseline' swSuffix]),...
            pDat.(['signifZBaseline_05' swSuffix]),...
            pDat.(['signifZBaseline_01' swSuffix])] = deal(rho,pval,sig05,sig01);
            % Use rhp/pval from Z-scored (with mean and std from) whole Trial period
            swRhoPval = pDat.(['rho_pval_static_Z_trial' swSuffix]);
            [rho,pval,sig05,sig01] = cellfun(@(x) deal(x(1),x(2),x(2)<=0.05,x(2)<=0.01),swRhoPval);            
            [pDat.(['rhoZTrial' swSuffix]),...
            pDat.(['pvalZTrial' swSuffix]),...
            pDat.(['signifZTrial_05' swSuffix]),...
            pDat.(['signifZTrial_01' swSuffix])] = deal(rho,pval,sig05,sig01);
        end
        
        %% xwaveforms...
        [pDat.xWaveMean,pDat.xWaveStd] = ...
            cellfun(@(w) deal(mean(cell2mat(w)),std(cell2mat(w))),pDat.xWaves, 'UniformOutput', false);
        [pDat.xWaveWidthMean,pDat.xWaveWidthStd] = ...
            cellfun(@(w) deal(mean(cell2mat(w)),std(cell2mat(w))),pDat.xWaveWidths, 'UniformOutput', false);
        %% ywaveforms...
        [pDat.yWaveMean,pDat.yWaveStd] = ...
            cellfun(@(w) deal(mean(cell2mat(w)),std(cell2mat(w))),pDat.yWaves, 'UniformOutput', false);
        [pDat.yWaveWidthMean,pDat.yWaveWidthStd] = ...
            cellfun(@(w) deal(mean(cell2mat(w)),std(cell2mat(w))),pDat.yWaveWidths, 'UniformOutput', false);
        
        outPairs(p).pDat = pDat;
    end
    
    spkCorrStatic.(areaPairField) = vertcat(outPairs.pDat); 
    spkCorrStatic.(areaPairField)(:,cols2Remove)=[];
    fprintf('Done %.3f sec.\n',toc)
end
size(spkCorrStatic);
[d,~,~]=fileparts(oFilename);
if ~exist(d,'dir')
    mkdir(d);
end


save(oFilename,'-v7.3','-struct','spkCorrStatic');
toc




