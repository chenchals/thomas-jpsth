
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




% choose unit no...
unitNo = 41;

unitName = num2str(unitNo,'Unit_%03d');
wavUnitTbl = load(fullfile(waveDataDir,[unitName '.mat']));
wavUnitTbl = wavUnitTbl.(unitName);
session = wavUnitTbl.session{1};

% All trials for that unit/session
spkUnitTrls = allSpks(ninfo.unitNum==unitNo).SAT';
% verify no of trials is correct for wav and spk
assert(size(wavUnitTbl.wavSearch{1},1)==numel(spkUnitTrls),'Error different no. of trials')

% Trial types for the session
trialTypes = sessionTrialTypes(strcmp(sessionTrialTypes.session,session),:);
% Event times
evntTimes = sessionEventTimes(strcmp(sessionEventTimes.session,session),:);
% for sel. trials get 'matching' waveforms for spks
selTrls = [10 12 128];
diffMax = 1; % 1 ms different
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
   alignedUnitTs = trlUnitTs;
   alignedWfTs = trlWfTs;
   
   % find indices of spikes that are closest match
   spkIdx = arrayfun(@(uTs) find(abs(alignedWfTs - uTs)<=diffMax,1),alignedUnitTs);
   plotWaveforms(trlWf(spkIdx,:))
   pause
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












