%% Variables for all units of SEF units paired with (SEF or FEF or SC) and compute SDFs
spkCorrDir = 'dataProcessed/analysis/11-18-2019/spkCorr';
datasetDir = 'dataProcessed/dataset';
sefPairsBaseFile = [spkCorrDir '/summary/sefPairs'];
oSdfPairsFile = fullfile(spkCorrDir,'summary/sefPairsUnitSdfs.mat');
spikeTimesFile = fullfile(datasetDir,'spikes_SAT.mat');
unitInfoStatsFile = fullfile(datasetDir,'dataNeurophys_SAT.mat');
trialTypesFile = fullfile(datasetDir,'TrialTypesDB.mat');
trialEventTimesFile = fullfile(datasetDir,'TrialEventTimesDB.mat');

%% load sefPairsBaseline and sefPairsPostSaccade files

blPairs = load([sefPairsBaseFile 'Baseline.mat'],'pairsByArea');
blPairs = blPairs.pairsByArea;
psPairs = load([sefPairsBaseFile 'PostSaccade.mat'],'pairsByArea');
psPairs = psPairs.pairsByArea;

%% Load other variables needed for computing SDFs
if ~exist('spikesSat','var')
    % Load data variable: spike times
    spikesSat = load(spikeTimesFile);
    spikesSat = spikesSat.spikesSAT;
    % Load unit Info for knowing which trials to be removed of any
    unitInfoAll = load(unitInfoStatsFile,'unitInfo');
    unitInfoAll = unitInfoAll.unitInfo;
    % Load data variable: TrialTypes
    sessionTrialTypes = load(trialTypesFile);
    sessionTrialTypes = sessionTrialTypes.TrialTypesDB;
    % Load data variable: TrialEventTimes
    sessionEventTimes = load(trialEventTimesFile);
    sessionEventTimes = sessionEventTimes.TrialEventTimesDB;
end


%% Compute SDFs for unique units for all pairs
if ~exist(oSdfPairsFile,'file')
    uniqUnits = unique([blPairs.unitNum;blPairs.pairedUnitNum;psPairs.unitNum;psPairs.pairedUnitNum]);
    % setup vars for SDF
    useAlignment.events = {'CueOn','SaccadePrimary','RewardTime'};
    useAlignment.timeWins = {[-600 400],[-200 600],[-100 400]};
    useAlignment.names = {'Visual','PostSaccade','PostReward'};
    useAlignment.sortEventNames = {'SaccadePrimary','SaccadePrimary','SaccadePrimary'};
    unitSdfsTbl = table();
    nUnits = numel(uniqUnits);
    str = '';
    for ii = 1:nUnits
        unitNum = uniqUnits(ii);
        fprintf(repmat('\b',1,length(str)+ 4));
        str = sprintf('Doing unit # %d [%d of %d] ...',unitNum,ii,nUnits);
        fprintf(str)
        unitInfo = unitInfoAll(unitNum,:);
        sess = unitInfo.sess{1};

        useUnit.unitNum = unitNum;
        useUnit.spkTimes = spikesSat{unitNum}';
        evntTimes = sessionEventTimes(strcmp(sessionEventTimes.session,sess),:);
        useTrials.trialTypes = sessionTrialTypes(strcmp(sessionTrialTypes.session,sess),:);
        useTrials.trialsToRemove = unitInfo.trRemSAT{1};
        useTrials.minTrialCount = 1;
        unitSdfsTbl = [unitSdfsTbl;getUnitSatSdf(useUnit,evntTimes,useTrials,useAlignment)]; %#ok<AGROW>
    end
    save(oSdfPairsFile,'-v7.3','unitSdfsTbl');
else
    load(oSdfPairsFile);
end
%% Normalize mean SDFs to maxFr = 1.0
fx_normMaxToOne = @(x) x./max(x); % max of sdf is 1.0 
% normalize unit SDFs
unitSdfsTbl.Visual_sdf = cellfun(@(x) fx_normMaxToOne(x(:,2))',unitSdfsTbl.Visual_sdfTsMeanStdSem,'UniformOutput',false);
unitSdfsTbl.PostSaccade_sdf = cellfun(@(x) fx_normMaxToOne(x(:,2))',unitSdfsTbl.PostSaccade_sdfTsMeanStdSem,'UniformOutput',false);
unitSdfsTbl.PostReward_sdf = cellfun(@(x) fx_normMaxToOne(x(:,2))',unitSdfsTbl.PostReward_sdfTsMeanStdSem,'UniformOutput',false);
% raw
unitSdfsTbl.Visual_raw_sdf = cellfun(@(x) x(:,2)',unitSdfsTbl.Visual_sdfTsMeanStdSem,'UniformOutput',false);
unitSdfsTbl.PostSaccade_raw_sdf = cellfun(@(x) x(:,2)',unitSdfsTbl.PostSaccade_sdfTsMeanStdSem,'UniformOutput',false);
unitSdfsTbl.PostReward_raw_sdf = cellfun(@(x) x(:,2)',unitSdfsTbl.PostReward_sdfTsMeanStdSem,'UniformOutput',false);
% ensure sdfs in Rows else convert
if ~size(cell2mat(unitSdfsTbl.Visual_sdf),1) == size(unitSdfsTbl,1)
    unitSdfsTbl.Visual_sdf = cellfun(@(x) x',unitSdfsTbl.Visual_sdf,'UniformOutput',false);
    unitSdfsTbl.PostSaccade_sdf = cellfun(@(x) x',unitSdfsTbl.PostSaccade_sdf,'UniformOutput',false);
    unitSdfsTbl.PostReward_sdf = cellfun(@(x) x',unitSdfsTbl.PostReward_sdf,'UniformOutput',false);
end

%% filter for pairArea (fieldnames)
pairedAreas = {{'SEF' 'SEF'};{'SEF' 'FEF'};{'SEF' 'SC'}};
% for SEF SEF pairs
% for SEF FEF pairs
% for SEF SC pairs

usePairs = [blPairs;psPairs];
for ii = 1:numel(pairedAreas)
    pairedArea = pairedAreas{ii};
    pairedAreaStr = char(join(pairedArea,'-'));
    area1 = pairedArea{1}; % SEF
    area2 = pairedArea{2}; % SEF or FEF or SC
    idxPair = strcmp(usePairs.X_area,area1) & strcmp(usePairs.Y_area,area2);
    
    units1 = unique(usePairs.unitNum(idxPair),'stable');
    units2 = unique(usePairs.pairedUnitNum(idxPair),'stable');
    
    units1Idx = find(ismember(unitSdfsTbl.unitNum,units1));
    units2Idx = find(ismember(unitSdfsTbl.unitNum,units2));
    %% 
    units1SdfsTbl = getSdfsTblForUnits(unitSdfsTbl,units1Idx,0,area1);
    units2SdfsTbl = getSdfsTblForUnits(unitSdfsTbl,units2Idx,0,area2);
    
    %%
    corrSpkSatSdfPlot(units1SdfsTbl,[],['AverageSdf_allPairedUnits_' pairedAreaStr '_a1_' area1 '.pdf']);
    corrSpkSatSdfPlot(units2SdfsTbl,[],['AverageSdf_allPairedUnits_' pairedAreaStr '_a2_' area2 '.pdf']);

    
    
end



%% Other functions
function [outSdfsTbl] = getSdfsTblForUnits(unitSdfsTbl,unitsIdx, unitNum, unitArea)
    % Use all pairs from baseline as well as postsaccade    
    fx_sdfTsMeanStdSem = @(t,x) {[...
        mean(cell2mat(t),1);...
        mean(cell2mat(x),1);...
        std(cell2mat(x),1);...
        std(cell2mat(x),1)./sqrt(size(x,1))...
        ]'};
    
    fx_allUnitsSdf = @(x) {cell2mat(x)}; % mean by rows


    temp = unitSdfsTbl(unitsIdx,:);
    [condGrpIdx,condGroupNames] = findgroups(temp.condition);

    outSdfsTbl = table();
    outSdfsTbl.condition = condGroupNames;
    % need sfds to be a matrix (nTs by 4), names _sdfTsMeanStdSem
    outSdfsTbl.Visual_sdfTsMeanStdSem = splitapply(@(t,x) fx_sdfTsMeanStdSem(t,x),temp.Visual_timeMs,temp.Visual_sdf,condGrpIdx);
    outSdfsTbl.PostSaccade_sdfTsMeanStdSem = splitapply(@(t,x) fx_sdfTsMeanStdSem(t,x),temp.PostSaccade_timeMs,temp.PostSaccade_sdf,condGrpIdx);
    outSdfsTbl.PostReward_sdfTsMeanStdSem = splitapply(@(t,x) fx_sdfTsMeanStdSem(t,x),temp.PostReward_timeMs,temp.PostReward_sdf,condGrpIdx);
    % add unitNum as 0 and area as specified in input
    outSdfsTbl.unitNum = repmat(unitNum,size(outSdfsTbl,1),1);
    outSdfsTbl.area = repmat({unitArea},size(outSdfsTbl,1),1);
    
    % all raw sdfs
    outSdfsTbl.Visual_raw_all = splitapply(fx_allUnitsSdf,temp.Visual_raw_sdf,condGrpIdx);
    outSdfsTbl.PostSaccade_raw_all = splitapply(fx_allUnitsSdf,temp.PostSaccade_raw_sdf,condGrpIdx);
    outSdfsTbl.PostReward_raw_all = splitapply(fx_allUnitsSdf,temp.PostReward_raw_sdf,condGrpIdx);    
    % all sdfs
    outSdfsTbl.Visual_all = splitapply(fx_allUnitsSdf,temp.Visual_sdf,condGrpIdx);
    outSdfsTbl.PostSaccade_all = splitapply(fx_allUnitsSdf,temp.PostSaccade_sdf,condGrpIdx);
    outSdfsTbl.PostReward_all = splitapply(fx_allUnitsSdf,temp.PostReward_sdf,condGrpIdx);

end
