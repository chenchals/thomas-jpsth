function [ results ] = analyzeBursts()
%%%%%%%%%%%%%%%%%%% User parameters Area %%%%%%%%%%%%%%%%%%%%%%%%
timeWin = [-1000 1500];
fileToLoad = '/Users/subravcr/Projects/lab-schall/schalllab-jpsth/data/spikeTimes_saccAligned_sess14.mat';
conditionsFile = '/Users/subravcr/Projects/lab-schall/schalllab-jpsth/data/ttx.mat';
session = 14;
condition = 'GO_Right';
% Load data
load(fileToLoad);
[~,datafile,ext] = fileparts(fileToLoad);
datafile = [datafile ext];
origSpkTimes = SpikeTimes.saccade;

% Unit Ids, channels, layers
unitIds = 1:size(origSpkTimes,2);
unitChannelIds = [1, 2, 3, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 9, 10, 10, 11, 12, 13, 14, 14, 14, 15, 17, 18, 19];
% Attached is the   ttx  file I explained to you in person.
% The session of interest is session 14, with 29 neurons.
% For the 29 neurons, the depths in channel units relative to the surface are:
% 1, 2, 3, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, ||||||| 9, 9, 9, 10, 10, 11, 12, 13, 14, 14, 14, 15, 17, 18, 19
%
% Conditions to select trials
load(conditionsFile);
[ ~, sessionConditionFile, ext] = fileparts(conditionsFile);
sessionConditionFile = [sessionConditionFile ext];
trials = ttx.(condition){session};

cellIdsTable = table();
cellIdsTable.unitIds = unitIds';
cellIdsTable.channelNo = unitChannelIds';
cellIdsTable.upperLayer = (unitChannelIds<=8)';
cellIdsTable.lowerLayer = (unitChannelIds>8)';

groups = {'units', 'channelUnits', 'layerUnits'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prune spkTimes . ?
spkTimesTrials = origSpkTimes(trials,:);
results = struct();
for g = 1:1 %numel(groups)
    
    clearvars prefix units2Use allPsth allPsthPsp psthBins allRasters rasterBins temp allBursts pltRows pltCols
    
    selGroup = groups{g};
    %units2Use = eval(selGroup);
    % use these Ids...
    switch selGroup
        case 'units'
            units2Use = arrayfun(@(x) {x},cellIdsTable.unitIds);
            titles = cellfun(@(x) num2str(x,'Unit Id #%02d'),units2Use,'UniformOutput',false);
            analysisType = 'Single Unit';
            pltRows = 5;
            pltCols = 6;
        case 'channelUnits'
            channels  = unique(cellIdsTable.channelNo);
            units2Use = arrayfun(@(x) find(cellIdsTable.channelNo==x),...
                unique(cellIdsTable.channelNo),'UniformOutput',false);
            titles = arrayfun(@(x) num2str(x,'Channel #%02d'),channels,'UniformOutput',false);
            analysisType = 'Pooled units per Channel';
            pltRows = 3;
            pltCols = 6;
        case 'layerUnits'
            units2Use = {
                % upper layer
                cell2mat(arrayfun(@(x) find(cellIdsTable.channelNo==x),...
                unique(cellIdsTable.channelNo(cellIdsTable.upperLayer)),'UniformOutput',false))';
                % lower layer
                cell2mat(arrayfun(@(x) find(cellIdsTable.channelNo==x),...
                unique(cellIdsTable.channelNo(cellIdsTable.lowerLayer)),'UniformOutput',false))';
                };
            titles = {'Upper layers' 'Lower layers'};
            analysisType = 'Pooled units per Layer';
            pltRows = 2;
            pltCols = 1;
            
    end
    
    % aggregate by channels
    spkTimes = SpikeUtils.groupSpikeTimes(spkTimesTrials, units2Use);
    
    % All Rasters
    temp = arrayfun(@(x) SpikeUtils.rasters(spkTimes(:,x),timeWin),...
        1:size(spkTimes,2),'UniformOutput',false);
    allRasters = cellfun(@(x) x.rasters,temp,'UniformOutput',false);
    rasterBins = temp{1}.rasterBins;
    clearvars temp
    
    % All PSTHs
    psthBinWidth = 1;
    temp = arrayfun(@(x) SpikeUtils.psth(spkTimes(:,x),psthBinWidth,timeWin,@nanmean),1:size(spkTimes,2),'UniformOutput',false);
    allPsth =  cellfun(@(x) x.psth,temp,'UniformOutput',false);
    psthBins = temp{1}.psthBins;
    allPsthPsp = cellfun(@(x) convn(x',SpikeUtils.pspKernel,'same')',allPsth,'UniformOutput',false);
    clearvars temp
    
    % Compute bursts
    plotProgress = 0;
    fprintf('Running burst detector...\n')
    %allBursts = cell(nTrials,numel(idsToUse));
    for c = 1:size(spkTimes,2)
        fprintf('%02d/%02d\n',c,size(spkTimes,2));
        x = spkTimes(:,c);
        allBursts(:,c) = BurstUtils.detectBursts(x,timeWin);
    end
    fprintf('done\n')
    plotResults(titles, allPsthPsp, psthBins, allRasters, rasterBins, allBursts, pltRows, pltCols);
    
    results(g).analysisTime = datetime;
    results(g).datafile = datafile;
    results(g).sessionConditionFile = sessionConditionFile;
    results(g).session = session;
    results(g).condition = condition;
    results(g).analysisPrefix = analysisType;
    results(g).titles = titles;
    results(g).unitGroups = units2Use;
    results(g).allPsthPsp = allPsthPsp;
    results(g).psthBins = psthBins;
    results(g).allRasters = allRasters;
    results(g).rasterBins = rasterBins;
    results(g).allBursts = allBursts;
    
    % further burst analysis
    fx_t = @(fn,cellBursts) cellfun(@(x) x.(fn),cellBursts,'UniformOutput',false);     
    for c = 1:size(allBursts,2)
        currBursts = allBursts(:,c);
        allBurstHists(:,c) = BurstUtils.psbh(fx_t('bobT',currBursts),fx_t('eobT',currBursts),timeWin);
    end
   
    results(g).allBurstHists = allBurstHists;
    
    
end
end


function plotResults(titles, psthAllGroups, psthTimes, rastersAllGroups, rasterTimes, burstsAllGroups, plotRows, plotCols)
%%%% Plot rasters %%%%
bobT = arrayfun(@(x) x{1}.bobT,burstsAllGroups,'UniformOutput',false);
eobT = arrayfun(@(x) x{1}.eobT,burstsAllGroups,'UniformOutput',false);
maxFrs = max(cellfun(@max,psthAllGroups));
maxFrRound = round(maxFrs/5)*5;

figure
for cellId = 1:numel(titles)
    currTitle = titles{cellId};
    subplot(plotRows,plotCols,cellId);
    PlotUtils.plotPsth(psthAllGroups{cellId}, psthTimes,maxFrRound);
    hold on
    PlotUtils.plotRastersAndBursts(rastersAllGroups{cellId},rasterTimes,bobT(:,cellId),eobT(:,cellId));
    hold off
    title(currTitle);
    xlabel('Saccade aligned');
    ylabel('Firing rate (Hz)');
    drawnow
end

end



