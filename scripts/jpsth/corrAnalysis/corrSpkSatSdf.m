% plot SDFs for different epochs of units used in spike count correlation
% for paired neurons

% Use static windows
fn = 'dataProcessed/analysis/11-18-2019/spkCorr/summary/spkCorrAllPairsStaticNew.mat';
%spkCorr = load(fn);
outPdfDir = 'dataProcessed/analysis/11-18-2019/spkCorr/summary/sdf';
if ~exist(outPdfDir,'dir')
    mkdir(outPdfDir);
end
outPdfFn = fullfile(outPdfDir,'satSdf_');



