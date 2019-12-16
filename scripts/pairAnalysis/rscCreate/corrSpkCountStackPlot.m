corrDirs = {
    'dataProcessed/analysis/spkCorr/spkCorr_SEF-SEF/mat' 
    'dataProcessed/analysis/spkCorr/spkCorr_SEF-FEF/mat'    
    'dataProcessed/analysis/spkCorr/spkCorr_SEF-SC/mat'     
    'dataProcessed/analysis/spkCorr/spkCorr_FEF-FEF/mat'    
    'dataProcessed/analysis/spkCorr/spkCorr_FEF-SC/mat'     
    'dataProcessed/analysis/spkCorr/spkCorr_SC-SC/mat'  
    };
plotDirs = regexprep(corrDirs,'/mat$','/pdf');
pairAreas = regexprep(corrDirs,{'.*_','/mat','-'},{'','','_'});

for p = 1:numel(plotDirs)
    if ~exist(plotDirs{p},'dir')
        mkdir(plotDirs{p})
    end
end


%%
corrInfoFields = {
    'condition'
    'alignedName'
    'alignedEvent'
    'alignedTimeWin'
    'alignTime'
    'rasterBins'
    'rho_pval_100ms'
    };
pairInfoFields = {
    'Pair_UID'
    'X_unitNum'
    'Y_unitNum'
    'X_area'
    'Y_area'
    'isOnSameChannel'
    };
%%
tic
oStruct = struct();
for d = 1:numel(corrDirs)
    %%
    pairArea = pairAreas{d};
    scFiles = dir([corrDirs{d},'/spkCorr*.mat']);
    scFiles = strcat({scFiles.folder}',filesep,{scFiles.name}');
    pdfDir = plotDirs{d};
    savePdfFlag = 1;
    %%
    %rpTemp = cell2table(cell(numel(scFiles),numel(pairInfoFields)+numel(corrInfoFields)));
    rpTemp = table();
    for f = 1:numel(scFiles)           
        %%
        scFile = scFiles{f};
        temp = load(scFile);
        % extract unitids for pair, rhoPval for 150ms moving win
        t = [temp.spkCorr(1,corrInfoFields) temp.cellPairInfo(1,pairInfoFields)]; 
        rpTemp = [rpTemp;[temp.spkCorr(1,corrInfoFields) temp.cellPairInfo(1,pairInfoFields)]];
     end
    %%
    oStruct.(pairArea) = rpTemp;
    
end
toc