% Extract coincidence and correlation histograms and other data from
% JPSTH-PAIR_nnnn.mat 
% Extracted data:

%%
jpsthDataBaseDir = 'dataProcessed/analysis/JPSTH';
outputDir = 'dataProcessed/analysis/JPSTH-extract1';

pairDirs = {
    'jpsth_SEF-FEF/mat'
    'jpsth_SEF-SC/mat'
    'jpsth_FEF-SC/mat'
    'jpsth_SEF-SEF/mat'
    'jpsth_FEF-FEF/mat'
    'jpsth_SC-SC/mat'
    };
outputFiles = cellfun(@(x) fileparts(x),pairDirs,'UniformOutput',false);

pairDirs = strcat(jpsthDataBaseDir,filesep,pairDirs);
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

xCorrWin = [-50 50];

% see doc pctRunOnAll
pctRunOnAll warning off;
for ii = 1:nDirs
    tic
    nRows = 0;
    pairDir = pairDirs{ii}; 
    fprintf('Doing pairDir %s\n',pairDir)
    outputFile = fullfile(outputDir,[lower(outputFiles{ii}) '.mat']);
    jpsthPairFiles = dir(fullfile(pairDir,'JPSTH-PAIR_*.mat'));
    jpsthPairFiles = strcat({jpsthPairFiles.folder}',filesep,{jpsthPairFiles.name}');
    cellPairInfo = table();
    outData = table();
    parfor jj = 1:numel(jpsthPairFiles)
        oD = table();
        pairFile = jpsthPairFiles{jj};
        fprintf('Processing file %s\n',pairFile)
        temp = load(pairFile,'cellPairInfo');
        cellPairInfo(jj,:) = temp.cellPairInfo;
        for kk = 1:numel(conditions)
            condition = conditions{kk};
            data = load(pairFile,condition);
            data = data.(condition);
            nRows = size(oD,1) + 1; 
            for ll = 1:numel(winChoices)
                oD.cellPairId{nRows} = temp.cellPairInfo.Pair_UID;
                oD.condition{nRows} = condition;
                rowName = winChoices(ll).name;
                winLims = winChoices(ll).winLims;
                oD.([rowName '_winLims']){nRows} = winLims;
                % 
                times = data{rowName,'xPsthBins'}{1}(:);
                startIdx = find(times==winLims(1));
                endIdx = find(times==winLims(2));                
                oD.([rowName '_coins']){nRows} = sum(data{rowName,'coincidenceHist'}{1}(startIdx:endIdx,2));
                oD.([rowName '_xSpikeCounts']){nRows} = sum(sum(data{rowName,'xSpikeCounts'}{1}(:,startIdx:endIdx)));
                oD.([rowName '_ySpikeCounts']){nRows} = sum(sum(data{rowName,'ySpikeCounts'}{1}(:,startIdx:endIdx)));
                %
                times = data{rowName,'xCorrHist'}{1}(:,1);               
                oD.([rowName '_xCorrWinLims']){nRows} = xCorrWin;
                oD.([rowName '_xCorrHist']){nRows} = sum(data{rowName,'xCorrHist'}{1}...
                                (times>=xCorrWin(1) & times<=xCorrWin(2),2)); %#ok<*PFBNS>
            end
        end
        outData = [outData;oD];
    end % end parfor 
    tempData.cellPairInfo = cellPairInfo;
    tempData.jpsthExtract = outData;
    save(outputFile,'-v7.3','-struct','tempData');
    fprintf('Done doing pairDir %s, took %d s\n',pairDir, round(toc,2))
end


