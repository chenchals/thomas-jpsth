
corrDirs = {
    'dataProcessed/analysis/spkCorr/spkCorr_SEF-SEF/mat' 
    'dataProcessed/analysis/spkCorr/spkCorr_SEF-FEF/mat'    
    'dataProcessed/analysis/spkCorr/spkCorr_SEF-SC/mat'     
    'dataProcessed/analysis/spkCorr/spkCorr_FEF-FEF/mat'    
    'dataProcessed/analysis/spkCorr/spkCorr_FEF-SC/mat'     
    'dataProcessed/analysis/spkCorr/spkCorr_SC-SC/mat'  
    };
%     'dataProcessed/analysis/spkCorr/spkCorr_SEF-NSEFN/mat'  
%     'dataProcessed/analysis/spkCorr/spkCorr_FEF-NSEFN/mat'  
%     'dataProcessed/analysis/spkCorr/spkCorr_SC-NSEFN/mat'   
%     'dataProcessed/analysis/spkCorr/spkCorr_NSEFN-NSEFN/mat'
%     };
plotDirs = regexprep(corrDirs,'/mat$','/pdf');

for p = 1:numel(plotDirs)
    if ~exist(plotDirs{p},'dir')
        mkdir(plotDirs{p})
    end
end

tic
for d = 1:numel(corrDirs)
    scFiles = dir([corrDirs{d},'/rscCorr*.mat']);
    scFiles = strcat({scFiles.folder}',filesep,{scFiles.name}');
    pdfDir = plotDirs{d};
    savePdfFlag = 1;
    for f = 1:numel(scFiles)      
        corrSpkCountPlot(scFiles{f},pdfDir,savePdfFlag);       
    end
end
toc
