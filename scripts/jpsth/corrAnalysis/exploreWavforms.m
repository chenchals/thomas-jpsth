
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
% inline function for aligning, windowing TS of trials
fx_alignTs = @(tsByTrl,alinTimeByTrl) arrayfun(@(idx) tsByTrl{idx} - alinTimeByTrl(idx),(1:numel(alinTimeByTrl))','UniformOutput',false);
fx_alindTsToWin = @(alindTsByTrl, tsWin) cellfun(@(x) x(x>=tsWin(1) & x<=tsWin(2)),alindTsByTrl,'UniformOutput',false);
% find matching indices for 2 arrrays of aligned ts
% max abs. diff is 1 ms
diffMax = 1; % 1 ms different
fx_matchedSpkIdx = @(alindUTs,alindWfTs) arrayfun(@(uTs) find(abs(alindWfTs - uTs)<=diffMax,1),alindUTs,'UniformOutput',false);



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







for tt = 1:numel(selTrls)
   trlNo = selTrls(tt);
   %trlStartTs = matData.TrialStart_(trlNo,1);
   % all waveforms for spikes in trial
   trlWf = wavUnitTbl.wavSearch{1}{trlNo};
   % timestamps got from waveforms
   trlWfTs =  wavUnitTbl.wavSearchTs{1}{trlNo};
   % timestamps from unit in mat file / ninfo
   trlUnitTs = spkUnitTrls{trlNo}';
   
   % Align and subset the spike times you want
   alindUnitTs = trlUnitTs;
   alignedWfTs = trlWfTs;
   
   % find indices of spikes that are closest match
   wfs{tt,1} = trlWf(fx_matchedSpkIdx(alindUnitTs,alignedWfTs),:);
   %plotWaveforms(trlWf(spkIdx,:))
   %pause
end


% Plot aggregated waveform for all selected trials and spikes
alignEvent = 'CueOn';

selTrls = find(trialTypes.FastErrorChoice{1});
alignedUnitTs = arrayfun(@(idx) spkUnitTrls{selTrls(idx)} -  





allWf = wavUnitTbl.wavSearch{1};
allWfTs = wavUnitTbl.wavSearchTs{1};
allUnitTs = [spkUnitTrls(selTrls)];
allWf4SelTrls = arrayfun(@(tn) {tn},selTrls);
allWfTs4SelTrls = arrayfun(@(tn) wavUnitTbl.wavSearch{1}{tn},selTrls);












