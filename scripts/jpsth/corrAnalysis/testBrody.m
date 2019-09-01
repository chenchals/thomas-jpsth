

%[thisCovariogram, sigHigh, sigLow, parts] = covariogramBrody(spike_1, spike_2, p1, p2, s1, s2, lag)
lag = 25;
spike_1 = double(Accurate.xRasters{1});
spike_2 = double(Accurate.yRasters{1});
p1 = Accurate.xPsth{1};
p2 = Accurate.yPsth{1};
s1 = Accurate.xPsthStd{1}.^2;
s2 = Accurate.yPsthStd{1}.^2;

vXAligned = Accurate.xCellSpikeTimes{1};
vYAligned = Accurate.yCellSpikeTimes{1};
vTimeWin = Accurate.alignedTimeWin{1};
binWidth = 5;

%[brody.thisCovariogram, brody.sigHigh, brody.sigLow, brody.parts] = covariogramBrody(spike_1, spike_2, p1, p2, s1, s2, lag);


jpsth = SpikeUtils.jpsth(vXAligned,vXAligned,vTimeWin,binWidth,50);

jer = SpikeUtils.jeromiahJpsth(vXAligned,vXAligned,vTimeWin,binWidth,50);

