wavDir = 'dataProcessed/dataset/waves2';

%% inline function for aligning, windowing TS of trials
fx_alignTs = @(tsByTrl,alinTimeByTrl) arrayfun(@(idx) tsByTrl{idx} - alinTimeByTrl(idx),(1:numel(alinTimeByTrl))','UniformOutput',false);
fx_alindTsToWin = @(alindTsByTrl, tsWin) cellfun(@(x) x(x>=tsWin(1) & x<=tsWin(2)),alindTsByTrl,'UniformOutput',false);
% find matching indices for 2 arrrays of aligned ts
% max abs. diff is 1 ms
diffMax = 1; % 1 ms different
fx_matchedSpkIdx = @(alindUTs,alindWfTs) arrayfun(@(uTs) find(abs(alindWfTs - uTs)<=diffMax,1),alindUTs,'UniformOutput',false);

%%
% X Unit
unitName = sprintf('Unit_%03d',cellPairInfo.X_unitNum);
allWavs = load(fullfile(wavDir,[unitName '.mat']));
wavTsX = allWavs.(unitName).wavSearchTs{1};
wavsX = allWavs.(unitName).wavSearch{1};
% Y Unit
unitName = sprintf('Unit_%03d',cellPairInfo.Y_unitNum);
allWavs = load(fullfile(wavDir,[unitName '.mat']));
wavTsY = allWavs.(unitName).wavSearchTs{1};
wavsY = allWavs.(unitName).wavSearch{1};


dat.xWaves = cell(12,1);
dat.yWaves = cell(12,1);

for jj = 1:size(dat,1)
selTrls = find(dat.trialNosByCondition{jj});
alignTime = dat.alignTime{jj};
alignWin = dat.alignedTimeWin{jj};
% align ts from waveform data
alindWfTsX = fx_alignTs(wavTsX(selTrls),alignTime);
alindWfTsXWin = fx_alindTsToWin(alindWfTsX,alignWin);
alindWfTsY = fx_alignTs(wavTsY(selTrls),alignTime);
alindWfTsYWin = fx_alindTsToWin(alindWfTsY,alignWin);
% aligned windowd unit ts
alindUnitTsXWin = dat.xCellSpikeTimes{jj};
alindUnitTsYWin = dat.yCellSpikeTimes{jj};
% find matching indices for extracting waveforms
spkIdxByTrlX = cellfun(@(aUts,aWts) fx_matchedSpkIdx(aUts,aWts),alindUnitTsXWin,alindWfTsXWin,'UniformOutput',false);
dat.xWaves{jj} = cellfun(@(wfByTrl,idxByTrl) wfByTrl(cell2mat(idxByTrl'),:),wavsX(selTrls),spkIdxByTrlX,'UniformOutput',false);
spkIdxByTrlY = cellfun(@(aUts,aWts) fx_matchedSpkIdx(aUts,aWts),alindUnitTsYWin,alindWfTsYWin,'UniformOutput',false);
dat.yWaves{jj} = cellfun(@(wfByTrl,idxByTrl) wfByTrl(cell2mat(idxByTrl'),:),wavsY(selTrls),spkIdxByTrlY,'UniformOutput',false);

end


