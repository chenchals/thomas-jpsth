%% Load spike corr matrix
spkCorrDir = 'dataProcessed/analysis/11-18-2019/spkCorr';
datasetDir = 'dataProcessed/dataset';
spkCorrFile = fullfile(spkCorrDir,'summary/spkCorrAllPairsStaticNew.mat');
spikeTimesFile = fullfile(datasetDir,'spikes_SAT.mat');
unitInfoStatsFile = fullfile(datasetDir,'dataNeurophys_SAT.mat');
trialTypesFile = fullfile(datasetDir,'TrialTypesDB.mat');
trialEventTimesFile = fullfile(datasetDir,'TrialEventTimesDB.mat');

%% filter for pairArea (fieldnames)
pairAreas = {'SEF_SEF','SEF_FEF','SEF_SC'}; % FEF_FEF, FEF_SC, SC_SC

rscPairArea = table();
for pa = pairAreas
    temp = load(spkCorrFile,char(pa));
    rscPairArea = [rscPairArea;temp.(char(pa))]; %#ok<*AGROW>
    clearvars temp
end



%% split by epoch

% for pairArea-epoch get spk corr distrib

% split distribution into thirds

% for each third - get units for one area and units for other area

% compute mean sdfs for all epochs

% z-score all mean sdfs with meanFr in a baseline window

% average all z-scored sdfs
