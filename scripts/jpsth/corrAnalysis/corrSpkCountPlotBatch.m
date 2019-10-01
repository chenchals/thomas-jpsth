
corrDirs = {
    'dataProcessed/analysis/SEF-PAPER/rSpkCounts/mat/CHOICE_ERR_PAIRS'
    'dataProcessed/analysis/SEF-PAPER/rSpkCounts/mat/TIMING_ERR_PAIRS'
    };
plotDirs = {
    'dataProcessed/analysis/SEF-PAPER/rSpkCounts/pdf/CHOICE_ERR_PAIRS'
    'dataProcessed/analysis/SEF-PAPER/rSpkCounts/pdf/TIMING_ERR_PAIRS'
    };
for p = 1:numel(plotDirs)
    if ~exist(plotDirs{p},'dir')
        mkdir(plotDirs{p})
    end
end

tic
for d = 1:2
    scFiles = dir([corrDirs{d},'/rscCorr*.mat']);
    scFiles = strcat({scFiles.folder}',filesep,{scFiles.name}');
    pdfDir = plotDirs{d};
    savePdfFlag = 1;
    for f = 1:numel(scFiles)      
        corrSpkCountPlot(scFiles{f},pdfDir,savePdfFlag);       
    end
end
toc
