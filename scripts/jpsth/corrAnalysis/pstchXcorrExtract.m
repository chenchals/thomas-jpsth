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

for ii = 1: 1 %nDirs
    nRows = 0;
    pairDir = pairDirs{ii}; 
    [~,outputFile] = fileparts(pairDir);
    outputFile = fullfile(outputDir,[outputFile '.mat']);
    jpsthPairFiles = dir(fullfile(pairDir,'JPSTH-PAIR_*.mat'));
    jpsthPairFiles = strcat({jpsthPairFiles.folder}',filesep,{jpsthPairFiles.name}');
    cellPairInfo = table();
    outData = table();
 
    for jj = 1:4 %numel(jpsthPairFiles)
        pairFile = jpsthPairFiles{jj};
        fprintf('Processing file %s\n',pairFile)
        temp = load(pairFile,'cellPairInfo');
        cellPairInfo(jj,:) = temp.cellPairInfo;
        for kk = 1:numel(conditions)
            condition = conditions{kk};
            data = load(pairFile,condition);
            data = data.(condition);
            nRows = size(outData,1) + 1;
            for ll = 1:numel(winChoices)
                outData.cellPairId{nRows} = temp.cellPairInfo.Pair_UID;
                outData.condition{nRows} = condition;
                rowName = winChoices(ll).name;
                winLims = winChoices(ll).winLims;
                outData.([rowName '_winLims']){nRows} = winLims;
                times = data{rowName,'xPsthBins'}{1}(:);
                startIdx = find(times==winLims(1));
                endIdx = find(times==winLims(2));
                
                outData.([rowName '_coins']){nRows} = sum(data{rowName,'coincidenceHist'}{1}(startIdx:endIdx,2));
            end
        end
        
    end
    
    
end


