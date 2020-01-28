% Create different types of SDF heatmaps 
% Data that needs to be loaded are for all pairs from
% summary/spkCorrAllPairsStaticRhoPval.mat 
%% Load static RSC matrix for all pairs of SEF units
    fn = 'dataProcessed/analysis/spkCorr/summary/spkCorrAllPairsStaticRhoPval.mat';
        % Load data
    temp = load(fn,'-regexp','SEF*');
    allSpkCorr = table();
    fns = fieldnames(temp);
    for jj = 1:numel(fns)
        allSpkCorr = [allSpkCorr;temp.(fns{jj})]; %#ok<*AGROW>
    end
    clearvars temp fns jj
%% Plot all units by significance FILTER: NONE for FR or SIGN of Corr
% Mean FR
[sdfsNoFilterMean,sdfsImgsNoFilterMean] = doNoFiltrMeanFr(allSpkCorr);
% Normalized FR
[sdfsNoFilter,sdfsImgsNoFilter] = doNoFiltrNormFr(allSpkCorr);
%
%% Plot all units by significance FILTER: BASELINE FIRING RATE for SEF and FEF cells
%    Set1: Baseline FR >= Median FR Baseline epoch; 
%    Set2: Baseline FR < Median FR Baseline epoch; 
[sdfsMedianFrBlGt,sdfsImgsMedianFrBlGt,...
    sdfsMedianFrBlLt,sdfsImgsMedianFrBlLt,allSpkCorr] = ...
    doMednFrBsLnEpochFiltrNorm(allSpkCorr);
%
%% Plot all units by significance FILTER: POST-SACCADE FIRING RATE for SEF and FEF cells
%    Set1: PostSaccade FR >= Median FR PostSacade epoch; 
%    Set2: PostSaccade FR < Median FR PostSacade epoch; 
[sdfsMedianFrPostSaccGt,sdfsImgsMedianFrPostSaccGt,...
    sdfsMedianFrPostSaccLt,sdfsImgsMedianFrPostSaccLt,allSpkCorr] = ...
    doMednFrPostSaccEpochFiltrNorm(allSpkCorr);


%% Plot all units by significance FILTER: POST-SACCADE FIRING RATE for SEF and FEF cells
%    Set1: PostSaccade FR >= Median FR PostSacade epoch; 
%    Set2: PostSaccade FR < Median FR PostSacade epoch; 







%% Extract into sub-functions for future plotting....

%% Plot all units by significance FILTER: POST-SACCADE FIRING RATE for SEF and FEF cells
%    Set1: PostSaccade FR >= Median FR PostSacade epoch; 
%    Set2: PostSaccade FR < Median FR PostSacade epoch; 
function [sdfsMedianFrPostSaccGt,sdfsImgsMedianFrPostSaccGt,...
          sdfsMedianFrPostSaccLt,sdfsImgsMedianFrPostSaccLt,allSpkCorr] = ...
          doMednFrPostSaccEpochFiltrNorm(allSpkCorr)
    parpoolSize = 0;
    if ~isempty(gcp('nocreate'))
        parpoolSize = Inf;
    end
    % get all units
    unitsList = unique([allSpkCorr.X_unitNum;allSpkCorr.Y_unitNum]);
    postSaccFrWin = [-100 200];
    sdfsDir = 'dataProcessed/dataset/satSdfs';
    unitsPostSaccFrTbl = table();
    % Compute mean PostSaccade firing rate for each unit
    parfor (un = 1:numel(unitsList),parpoolSize)
        fr = table();
        unitNum = unitsList(un);
        fn = fullfile(sdfsDir,num2str(unitNum,'Unit_%03d.mat'));
        temp = load(fn,'sdfs');
        ts = temp.sdfs.PostSaccade_timeMs{1};
        idxX = find(ts>= postSaccFrWin(1) & ts <= postSaccFrWin(2)); %#ok<*PFBNS>
        temp = cat(1,temp.sdfs.PostSaccade_sdfByTrial{:});
        temp = temp(:,idxX);
        fr.unitNum = unitNum;
        fr.meanPostSaccFr = mean(mean(temp,2));
        unitsPostSaccFrTbl = [unitsPostSaccFrTbl;fr];
    end

    % add units' mean PostSaccade FR to the allSpkCorr X_units and Y_units
    % do for X_unitNums
    allSpkCorr=innerjoin(allSpkCorr,unitsPostSaccFrTbl,'LeftKeys','X_unitNum','RightKeys','unitNum',...
        'RightVariables','meanPostSaccFr');
    allSpkCorr.X_meanPostSaccFr = allSpkCorr.meanPostSaccFr;
    allSpkCorr.meanPostSaccFr = [];
    % do for Y_unitNums
    allSpkCorr=innerjoin(allSpkCorr,unitsPostSaccFrTbl,'LeftKeys','Y_unitNum','RightKeys','unitNum',...
        'RightVariables','meanPostSaccFr');
    allSpkCorr.Y_meanPostSaccFr = allSpkCorr.meanPostSaccFr;
    allSpkCorr.meanPostSaccFr = [];
    % Find median Baseline FR for SEF, FEF, SC
    unitAreas = {'SEF','FEF','SC'};
    medianPostSaccFr = struct();
    for aa = 1:numel(unitAreas)
        unitArea = unitAreas{aa};
        idxX = ismember(allSpkCorr.X_area,unitArea);
        idxY = ismember(allSpkCorr.Y_area,unitArea);
        medianPostSaccFr.(unitArea) = median([unique(allSpkCorr.X_meanPostSaccFr(idxX),'stable');...
            unique(allSpkCorr.Y_meanPostSaccFr(idxY),'stable')]);
    end
    % Filter allSpkCorr by median PostSaccade FR for each type of pairs
    % SEF-FEF pairs filter for both SEF and FEF
    idxGt_Sef_Fef = allSpkCorr.X_meanPostSaccFr >= medianPostSaccFr.SEF & allSpkCorr.Y_meanPostSaccFr >= medianPostSaccFr.FEF ...
        & ismember(allSpkCorr.X_area,'SEF') & ismember(allSpkCorr.Y_area,'FEF');
    idxLt_Sef_Fef = allSpkCorr.X_meanPostSaccFr < medianPostSaccFr.SEF & allSpkCorr.Y_meanPostSaccFr < medianPostSaccFr.FEF ...
        & ismember(allSpkCorr.X_area,'SEF') & ismember(allSpkCorr.Y_area,'FEF');
    % SEF-SC pairs filter only for SEF and NOT SC
    idxGt_Sef_Sc = allSpkCorr.X_meanPostSaccFr >= medianPostSaccFr.SEF...
        & ismember(allSpkCorr.X_area,'SEF') & ismember(allSpkCorr.Y_area,'SC');
    idxLt_Sef_Sc = allSpkCorr.X_meanPostSaccFr < medianPostSaccFr.SEF...
        & ismember(allSpkCorr.X_area,'SEF') & ismember(allSpkCorr.Y_area,'SC');
    % SEF-SEF pairs filter for both X_SEF and Y_SEF
    idxGt_Sef_Sef = allSpkCorr.X_meanPostSaccFr >= medianPostSaccFr.SEF & allSpkCorr.Y_meanPostSaccFr >= medianPostSaccFr.SEF ...
        & ismember(allSpkCorr.X_area,'SEF') & ismember(allSpkCorr.Y_area,'SEF');
    idxLt_Sef_Sef = allSpkCorr.X_meanPostSaccFr < medianPostSaccFr.SEF & allSpkCorr.Y_meanPostSaccFr < medianPostSaccFr.SEF ...
        & ismember(allSpkCorr.X_area,'SEF') & ismember(allSpkCorr.Y_area,'SEF');
    medianPostSaccFrGtTbl = [allSpkCorr(idxGt_Sef_Fef,:);allSpkCorr(idxGt_Sef_Sc,:);allSpkCorr(idxGt_Sef_Sef,:)];
    medianPostSaccFrLtTbl = [allSpkCorr(idxLt_Sef_Fef,:);allSpkCorr(idxLt_Sef_Sc,:);allSpkCorr(idxLt_Sef_Sef,:)];
    % process heatmaps for median gt table
    unitsMedianFrGtTbl = categorizeUnitsByRscSignif(medianPostSaccFrGtTbl);
    useNormalized = 1; subDirName = 'medianPostSaccFrGtNorm';
    [sdfsMedianFrPostSaccGt,sdfsImgsMedianFrPostSaccGt] = plotHeatmapPanels(unitsMedianFrGtTbl,useNormalized,subDirName);
    % process heatmaps for median lt table
    unitsMedianFrLtTbl = categorizeUnitsByRscSignif(medianPostSaccFrLtTbl);
    useNormalized = 1; subDirName = 'medianPostSaccFrLtNorm';
    [sdfsMedianFrPostSaccLt,sdfsImgsMedianFrPostSaccLt] = plotHeatmapPanels(unitsMedianFrLtTbl,useNormalized,subDirName);

end

%% Plot all units by significance FILTER: BASELINE FIRING RATE for SEF and FEF cells
%    Set1: Baseline FR >= Median FR Baseline epoch; 
%    Set2: Baseline FR < Median FR Baseline epoch; 
function [sdfsMedianFrBlGt,sdfsImgsMedianFrBlGt,sdfsMedianFrBlLt,sdfsImgsMedianFrBlLt,allSpkCorr] = ...
        doMednFrBsLnEpochFiltrNorm(allSpkCorr)

    parpoolSize = 0;
    if ~isempty(gcp('nocreate'))
        parpoolSize = Inf;
    end
    % get all units
    unitsList = unique([allSpkCorr.X_unitNum;allSpkCorr.Y_unitNum]);
    blFrWin = [-600 0];
    sdfsDir = 'dataProcessed/dataset/satSdfs';
    unitsFrTbl = table();
    % Compute mean baseline firing rate for each unit
    parfor (un = 1:numel(unitsList),parpoolSize)
        fr = table();
        unitNum = unitsList(un);
        fn = fullfile(sdfsDir,num2str(unitNum,'Unit_%03d.mat'));
        temp = load(fn,'sdfs');
        ts = temp.sdfs.Visual_timeMs{1};
        idxX = find(ts>= blFrWin(1) & ts <= blFrWin(2)); %#ok<*PFBNS>
        temp = cat(1,temp.sdfs.Visual_sdfByTrial{:});
        temp = temp(:,idxX); %#ok<*FNDSB>
        fr.unitNum = unitNum;
        fr.meanBlFr = mean(mean(temp,2));
        unitsFrTbl = [unitsFrTbl;fr];
    end

    % add units' mean Baseline FR to the allSpkCorr X_units and Y_units
    % do for X_unitNums
    allSpkCorr=innerjoin(allSpkCorr,unitsFrTbl,'LeftKeys','X_unitNum','RightKeys','unitNum',...
        'RightVariables','meanBlFr');
    allSpkCorr.X_meanBlFr = tempX.meanBlFr;
    allSpkCorr.meanBlFr = [];
    % do for Y_unitNums
    allSpkCorr=innerjoin(allSpkCorr,unitsFrTbl,'LeftKeys','Y_unitNum','RightKeys','unitNum',...
        'RightVariables','meanBlFr');
    allSpkCorr.Y_meanBlFr = allSpkCorr.meanBlFr;
    allSpkCorr.meanBlFr = [];
    % Find median Baseline FR for SEF, FEF, SC
    unitAreas = {'SEF','FEF','SC'};
    medianBlFr = struct();
    for aa = 1:numel(unitAreas)
        unitArea = unitAreas{aa};
        idxX = ismember(allSpkCorr.X_area,unitArea);
        idxY = ismember(allSpkCorr.Y_area,unitArea);
        %unitNums = unique(allSpkCorr.X_unitNum(idx),'stable');
        %zz= arrayfun(@(x) unique(allSpkCorr.X_meanBlFr(allSpkCorr.X_unitNum==x)),unitNums);
        medianBlFr.(unitArea) = median([unique(allSpkCorr.X_meanBlFr(idxX),'stable');...
            unique(allSpkCorr.Y_meanBlFr(idxY),'stable')]);
    end

    % Filter allSpkCorr by median balseline FR for each type of pairs
    % SEF-FEF pairs filter for both SEF and FEF
    idxGt_Sef_Fef = allSpkCorr.X_meanBlFr >= medianBlFr.SEF & allSpkCorr.Y_meanBlFr >= medianBlFr.FEF ...
        & ismember(allSpkCorr.X_area,'SEF') & ismember(allSpkCorr.Y_area,'FEF');
    idxLt_Sef_Fef = allSpkCorr.X_meanBlFr < medianBlFr.SEF & allSpkCorr.Y_meanBlFr < medianBlFr.FEF ...
        & ismember(allSpkCorr.X_area,'SEF') & ismember(allSpkCorr.Y_area,'FEF');
    % SEF-SC pairs filter only for SEF and NOT SC
    idxGt_Sef_Sc = allSpkCorr.X_meanBlFr >= medianBlFr.SEF...
        & ismember(allSpkCorr.X_area,'SEF') & ismember(allSpkCorr.Y_area,'SC');
    idxLt_Sef_Sc = allSpkCorr.X_meanBlFr < medianBlFr.SEF...
        & ismember(allSpkCorr.X_area,'SEF') & ismember(allSpkCorr.Y_area,'SC');
    % SEF-SEF pairs filter for both X_SEF and Y_SEF
    idxGt_Sef_Sef = allSpkCorr.X_meanBlFr >= medianBlFr.SEF & allSpkCorr.Y_meanBlFr >= medianBlFr.SEF ...
        & ismember(allSpkCorr.X_area,'SEF') & ismember(allSpkCorr.Y_area,'SEF');
    idxLt_Sef_Sef = allSpkCorr.X_meanBlFr < medianBlFr.SEF & allSpkCorr.Y_meanBlFr < medianBlFr.SEF ...
        & ismember(allSpkCorr.X_area,'SEF') & ismember(allSpkCorr.Y_area,'SEF');
    medianFrGtTbl = [allSpkCorr(idxGt_Sef_Fef,:);allSpkCorr(idxGt_Sef_Sc,:);allSpkCorr(idxGt_Sef_Sef,:)];
    medianFrLtTbl = [allSpkCorr(idxLt_Sef_Fef,:);allSpkCorr(idxLt_Sef_Sc,:);allSpkCorr(idxLt_Sef_Sef,:)];
    % process heatmaps for median gt table
    unitsMedianFrGtTbl = categorizeUnitsByRscSignif(medianFrGtTbl);
    useNormalized = 1; subDirName = 'medianFrBaselineGtNorm';
    [sdfsMedianFrBlGt,sdfsImgsMedianFrBlGt] = plotHeatmapPanels(unitsMedianFrGtTbl,useNormalized,subDirName);
    % process heatmaps for median lt table
    unitsMedianFrLtTbl = categorizeUnitsByRscSignif(medianFrLtTbl);
    useNormalized = 1; subDirName = 'medianFrBaselineLtNorm';
    [sdfsMedianFrBlLt,sdfsImgsMedianFrBlLt] = plotHeatmapPanels(unitsMedianFrLtTbl,useNormalized,subDirName);

end

%% Plot all units by significance FILTER: NONE for FR or SIGN of Corr
function [sdfsNoFilterMean,sdfsImgsNoFilterMean] = doNoFiltrMeanFr(allSpkCorr)
unitsNoFilterTbl = categorizeUnitsByRscSignif(allSpkCorr);
useNormalized = 0; subDirName = 'noFilterMean';
[sdfsNoFilterMean,sdfsImgsNoFilterMean] = plotHeatmapPanels(unitsNoFilterTbl,useNormalized,subDirName);
end

function [sdfsNoFilter,sdfsImgsNoFilter] = doNoFiltrNormFr(allSpkCorr)
unitsNoFilterTbl = categorizeUnitsByRscSignif(allSpkCorr);
useNormalized = 1; subDirName = 'noFilterNorm';
[sdfsNoFilter,sdfsImgsNoFilter] = plotHeatmapPanels(unitsNoFilterTbl,useNormalized,subDirName);
end

