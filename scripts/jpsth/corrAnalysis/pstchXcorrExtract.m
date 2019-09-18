% Extract coincidence and correlation histograms and other data from
% JPSTH-PAIR_nnnn.mat 
% Extracted data:

%%
jpsthDataBaseDir = 'dataProcessed/analysis/JPSTH';

pairDirs = {
    'jpsth_SEF-FEF/mat'
    'jpsth_SEF-SC/mat'
    'jpsth_FEF-SC/mat'
    'jpsth_SEF-SEF/mat'
    'jpsth_FEF-FEF/mat'
    'jpsth_SC-SC/mat'
    };

pairDirs = strcat(jpsthDataBaseDir,filesep,pairDirs);

outputDir = 'dataProcessed/analysis/JPSTH-extract1';
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

nDirs = numel(pairDirs);

conditions = {
    'Accurate'
    'AccurateCorrect'
    'AccurateErrorChoice'
    'AccurateErrorTiming'
    'Fast'
    'FastCorrect'
    'FastErrorChoice'
    'FastErrorTiming'
    };
winChoices(1).name = 'Baseline';
winChoices(1).winLims = [-500 -100];
winChoices(2).name = 'Visual';
winChoices(2).winLims = [50 200];
winChoices(3).name = 'PostSaccade';
winChoices(3).winLims = [50 300];


for ii = 1:nDirs
    pairDir = pairDirs{ii}; 
    [~,outputFile] = fileparts(pairDir);
    outputFile = fullfile(outputDir,[outputFile '.mat']);
    jpsthPairFiles = dir(fullfile(pairDir,'JPSTH-PAIR_*.mat'));
    jpsthPairFiles = strcat({jpsthPairFiles.folder}',filesep,{jpsthPairFiles.name}');
    cellPairInfo = table();
    outData = table();
 
    parfor jj = 1:4 %numel(jpsthPairFiles)
        pairFile = jpsthPairFiles{jj};
        fprintf('Processing file %s\n',pairFile)
        temp = load(pairFile,'cellPairInfo');
        cellPairInfo(jj,:) = temp.cellPairInfo;
        nRows = size(outData,1) + 1;
        for kk = 1:numel(conditions)
            nRows = nRows + 1;
            condition = conditions{kk};
            data = load(pairFile,condition);
            data = data.(condition);
            outData(nRows).condition = condition;
            for ll = 1:numel(winChoices)
                rowName = winChoices(ll).name;
                winLims = winChoices(ll).winLims;
                outData(jj).([rowName '_winLims']) = winLims;
                times = data{rowName,'xPsthBins'}{1}(:);
                startIdx = find(times==winLims(1));
                endIdx = find(times==winLims(2));
                
                outData(jj).([rowName '_coins']) = sum(data{rowName,'coincidenceHist'}{1}(startIdx:endIdx,2));
            end
        end
        
    end
    
    
end


