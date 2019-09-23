
% matfile 'teba/data' dir mapping location
matDataDir = 'T:'; 
% load all waveforms for SAT
%allUnits = load('E:/waves_SAT_DaEu1.mat');
%allUnits = load('dataProcessed/dataset/waves_SAT_DaEu1.mat');

% list of units:
unitList = fieldnames(allUnits);
% remove all names that do not start with Unit_
% cheat...
unitList = unitList(1:end-2);

% load all SAT spikes
allSpks = load('dataProcessed\dataset\spikes_SAT.mat');
allSpks = allSpks.spikes;

% load ninfo for getting the unit number / unit name...
ninfo = load('dataProcessed\dataset\ninfo_nstats_SAT.mat','ninfo');
ninfo = struct2table(ninfo.ninfo);
% choose unit no...
unitNo = 41;
wavUnitTbl = allUnits.(num2str(unitNo,'Unit_%03d'));
spkUnitTrls = allSpks(ninfo.unitNum==unitNo).SAT';
% verify no of trials is correct for wav and spk
assert(size(wavUnitTbl.wavSearch{1},1)==numel(spkUnitTrls),'Error different no. of trials')

% get the TrialStart_ time from the mat file...
matData = load(fullfile(matDataDir,wavUnitTbl.matFile{1}),'TrialStart_');



% for sel. trials get 'matching' waveforms for spks
selTrls = [10 12 128];
diffMax = 1; % 1 ms different
for tt = 1:numel(selTrls)
   trlNo = selTrls(tt);
   trlStartTs = matData.TrialStart_(trlNo,1);
   % all waveforms for spikes in trial
   trlWf = wavUnitTbl.wavSearch{1}{trlNo};
   % timestamps got from waveforms
   trlWfTs =  round((wavUnitTbl.wavSearchTs{1}{trlNo}-trlStartTs));
   % timestamps from unit in mat file / ninfo
   trlUnitTs = spkUnitTrls{trlNo}';
   % find indices of spikes that are closest match
   spkIdx = arrayfun(@(uTs) find(abs(trlWfTs - uTs)<=diffMax,1),trlUnitTs);
   plotWaveforms(trlWf(spkIdx,:))
   
   
   plot(1:32,trlWf(1:end,:))
end












