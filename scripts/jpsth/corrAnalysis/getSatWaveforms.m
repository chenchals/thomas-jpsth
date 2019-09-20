jpsthPairs = load('dataProcessed/dataset/JPSTH_PAIRS_CellInfoDB.mat');
jpsthPairs = jpsthPairs.JpsthPairCellInfoDB;
ninfo = load('dataProcessed/dataset/ninfo_nstats_SAT.mat');
ninfo = struct2table(ninfo.ninfo);
%Choice Errors (n=24): find( [ninfo.errGrade] >= 2)
eChoice = ninfo(ninfo.errGrade >= 2,:);
jpsthPairs.tempUnitNum = jpsthPairs.X_unitNum;
A = innerjoin(eChoice,jpsthPairs,'LeftKeys','unitNum','RightKeys',{'tempUnitNum'});
jpsthPairs.tempUnitNum = jpsthPairs.Y_unitNum;
B = innerjoin(eChoice,jpsthPairs,'LeftKeys','unitNum','RightKeys',{'tempUnitNum'});
eChoicePairs = sortrows([A;B],'unitNum');
clearvars A B

%Timing Errors (n = 25): find( abs([ninfo.rewGrade]) >= 2)
eTiming = ninfo(abs(ninfo.rewGrade) >= 2,:);
jpsthPairs.tempUnitNum = jpsthPairs.X_unitNum;
A = innerjoin(eTiming,jpsthPairs,'LeftKeys','unitNum','RightKeys',{'tempUnitNum'});
jpsthPairs.tempUnitNum = jpsthPairs.Y_unitNum;
B = innerjoin(eTiming,jpsthPairs,'LeftKeys','unitNum','RightKeys',{'tempUnitNum'});
eTtmingPairs = sortrows([A;B],'unitNum');

clearvars A B jpsthPairs


%% Error Choice Pair 0189 - 13a, 13b
pairInfo = eChoicePairs(strcmpi(eChoicePairs.Pair_UID,'PAIR_0189'),:);
matFile = pairInfo.matDatafile{1}; % maybe duplicates are present
% need ot translate or do manually for all sessions
plxFile = 'data/Euler/SAT/Plexon/Sorted/E20130911001-RH.plx';
plxData = readPLXFileC(plxFile,'all');
chData = plxData.SpikeChannels([plxData.SpikeChannels.Channel]==13);
spikeData.matFile = matFile;
spikeData.plexonFile = 'data/Euler/SAT/Plexon/Sorted/E20130911001-RH.plx';
spikeData.pairID = pairInfo.Pair_UID;
spikeData.channel = 13;
spikeData.clustNo = unique(chData.Units);
spikeData.spkCounts = arrayfun(@(x) sum(chData.Units==x),spikeData.clustNo,'UniformOutput',false);
spikeData.spkTimes = arrayfun(@(x) double(chData.Timestamps(chData.Units==x)),spikeData.clustNo,'UniformOutput',false);
spikeData.waves = arrayfun(@(x) double(chData.Waves(:,chData.Units==x)),spikeData.clustNo,'UniformOutput',false);
spikeData.mean = cellfun(@(x) mean(x,2),spikeData.waves,'UniformOutput',false);
spikeData.std = cellfun(@(x) std(x,[],2),spikeData.waves,'UniformOutput',false);
plotWfMeanVar(spikeData);


%%





function plotWfMeanVar(spikeData)
    colors = {'b','r','c'};
    for ii = 1:numel(spikeData.mean)
        spkCounts = spikeData.spkCounts{ii};
        clustNum = spikeData.clustNo(ii);
        clustLetter = char(clustNum+97);
        labels{2*ii-1} = '';
        labels{2*ii} = sprintf('Unit %d%s n=%d',spikeData.channel,clustLetter,spkCounts);
        y = (spikeData.mean{ii})';
        x = (1:32);
        yStd = (spikeData.std{ii})';
        fill([x fliplr(x)], [y+yStd fliplr(y-yStd)],colors{ii},'FaceAlpha',0.5,'LineStyle','none')
        hold on
        plot(x,y,colors{ii},'LineWidth',1.5)
    end
    legend(labels,'Location','northwest','Box','off')
    hold off
end

function plotWf(wavformMat)
% Each column is a waveform
  nPoints = size(wavformMat,1);
  nWaves = size(wavformMat,2);
  temp = [wavformMat;nan(1,nWaves)];
  temp = temp(:);
  x = repmat([1:nPoints NaN]',1,nWaves);
  x = x(:);
  plot(x,temp)
end

