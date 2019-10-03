
waveDataDir = 'dataProcessed/dataset/waves';
% load all SAT spikes
allSpks = load('dataProcessed/dataset/spikes_SAT.mat');
allSpks = allSpks.spikes;

% load ninfo for getting the unit number / unit name...
ninfo = load('dataProcessed/dataset/ninfo_nstats_SAT.mat','ninfo');
ninfo = struct2table(ninfo.ninfo);

% get session event times
sessionEventTimes = load('dataProcessed/dataset/TrialEventTimesDB.mat');
sessionEventTimes = sessionEventTimes.TrialEventTimesDB;

% get Conditions for sessions (Accurate, Fast, ...)
sessionTrialTypes = load('dataProcessed/dataset/TrialTypesDB.mat');
sessionTrialTypes = sessionTrialTypes.TrialTypesDB;

%% inline function for aligning, windowing TS of trials
fx_alignTs = @(tsByTrl,alinTimeByTrl) arrayfun(@(idx) tsByTrl{idx} - alinTimeByTrl(idx),(1:numel(alinTimeByTrl))','UniformOutput',false);
fx_alindTsToWin = @(alindTsByTrl, tsWin) cellfun(@(x) x(x>=tsWin(1) & x<=tsWin(2)),alindTsByTrl,'UniformOutput',false);
% find matching indices for 2 arrrays of aligned ts
% max abs. diff is 1 ms
diffMax = 1; % 1 ms different
fx_matchedSpkIdx = @(alindUTs,alindWfTs) arrayfun(@(uTs) find(abs(alindWfTs - uTs)<=diffMax,1),alindUTs,'UniformOutput',false);


%% for a given unit, condition, alignment, align-window plot selected waveforms
% choose unit no...
unitNo = 41;
% trialType
selTrialType = 'FastErrorChoice';
% note CueOn == TrialStart
alignedEvent = 'CueOn';
alignWin = [-200 400];

% load unit
unitName = num2str(unitNo,'Unit_%03d');
wavUnitTbl = load(fullfile(waveDataDir,[unitName '.mat']));
wavUnitTbl = wavUnitTbl.(unitName);
session = wavUnitTbl.session{1};

% Event times
evntTimes = sessionEventTimes(strcmp(sessionEventTimes.session,session),:);
alignTime = evntTimes.CueOn{1};
if ~strcmp(alignedEvent,'CueOn')
    alignTime = alignTime + evntTimes.(alignedEvent){1}(:);
end

% load trials for that unit/session
spkUnitTrls = allSpks(ninfo.unitNum==unitNo).SAT';
% verify no of trials is correct for wav and spk
assert(size(wavUnitTbl.wavSearch{1},1)==numel(spkUnitTrls),'Error different no. of trials')

% Selected Trial types for the session
trialTypes = sessionTrialTypes(strcmp(sessionTrialTypes.session,session),:);
% trial selection, and alignment
selTrls = find(trialTypes.(selTrialType){1});
% align and select ts in window after alignment
alindUnitTs = fx_alignTs(spkUnitTrls(selTrls),alignTime(selTrls));
alindUnitTsWin = fx_alindTsToWin(alindUnitTs,alignWin);
% align ts from waveform data
alindWfTs = fx_alignTs(wavUnitTbl.wavSearchTs{1}(selTrls),alignTime(selTrls));
alindWfTsWin = fx_alindTsToWin(alindWfTs,alignWin);
% find matching indices for extracting waveforms
spkIdxByTrl = cellfun(@(aUts,aWts) fx_matchedSpkIdx(aUts,aWts),alindUnitTsWin,alindWfTsWin,'UniformOutput',false);
% all waveforms for spikes in selected trials
allWfBytrl = wavUnitTbl.wavSearch{1}(selTrls);
alindWfWin = cellfun(@(wfByTrl,idxByTrl) wfByTrl(cell2mat(idxByTrl'),:),allWfBytrl,spkIdxByTrl,'UniformOutput',false);

plotWaveforms(cell2mat(alindWfWin));











