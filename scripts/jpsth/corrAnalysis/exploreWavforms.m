
waveDataDir = 'dataProcessed/dataset/waves';
% load all SAT spikes
allSpks = load('dataProcessed/dataset/spikes_SAT.mat');
allSpks = allSpks.spikes;

% load ninfo for getting the unit number / unit name...
ninfo = load('dataProcessed/dataset/ninfo_nstats_SAT.mat','ninfo');
ninfo = struct2table(ninfo.ninfo);
% choose unit no...
unitNo = 41;

unitName = num2str(unitNo,'Unit_%03d');
wavUnitTbl = load(fullfile(waveDataDir,[unitName '.mat']));
wavUnitTbl = wavUnitTbl.(unitName);
spkUnitTrls = allSpks(ninfo.unitNum==unitNo).SAT';
% verify no of trials is correct for wav and spk
assert(size(wavUnitTbl.wavSearch{1},1)==numel(spkUnitTrls),'Error different no. of trials')

% for sel. trials get 'matching' waveforms for spks
selTrls = [10 12 128];
diffMax = 1; % 1 ms different
for tt = 1:numel(selTrls)
   trlNo = selTrls(tt);
   trlStartTs = matData.TrialStart_(trlNo,1);
   % all waveforms for spikes in trial
   trlWf = wavUnitTbl.wavSearch{1}{trlNo};
   % timestamps got from waveforms
   trlWfTs =  wavUnitTbl.wavSearchTs{1}{trlNo};
   % timestamps from unit in mat file / ninfo
   trlUnitTs = spkUnitTrls{trlNo}';
   % find indices of spikes that are closest match
   spkIdx = arrayfun(@(uTs) find(abs(trlWfTs - uTs)<=diffMax,1),trlUnitTs);
   plotWaveforms(trlWf(spkIdx,:))
end












