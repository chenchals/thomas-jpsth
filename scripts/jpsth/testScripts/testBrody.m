

%[thisCovariogram, sigHigh, sigLow, parts] = covariogramBrody(spike_1, spike_2, p1, p2, s1, s2, lag)
% jer_default_lag = 25;
% spike_1 = double(Accurate.xRasters{1});
% spike_2 = double(Accurate.yRasters{1});
% p1 = Accurate.xPsth{1};
% p2 = Accurate.yPsth{1};
% s1 = Accurate.xPsthStd{1}.^2;
% s2 = Accurate.yPsthStd{1}.^2;
%%
vXAligned = AccurateErrorTiming.xCellSpikeTimes{3};
vYAligned = AccurateErrorTiming.yCellSpikeTimes{3};
vTimeWin = AccurateErrorTiming.alignedTimeWin{3};

%[brody.thisCovariogram, brody.sigHigh, brody.sigLow, brody.parts] = covariogramBrody(spike_1, spike_2, p1, p2, s1, s2, lag);
binWidth = 1;
coinsBinWidth = 50;
jer_default_lag = 50;

jpsth = SpikeUtils.jpsth(vXAligned,vYAligned,vTimeWin,binWidth,coinsBinWidth);

jer = SpikeUtils.jeromiahJpsth(vXAligned,vYAligned,vTimeWin,binWidth,coinsBinWidth);
%%
figure
subplot(2,4,1)
imagesc(flipud(jpsth.normalizedJpsth))
title('Jer-Normalized JPSTH')

subplot(2,4,2)
area(jpsth.coincidenceHist(:,1),smoothdata(jpsth.coincidenceHist(:,2),'gaussian' ,25,'omitnan'),'EdgeColor','none')

title(sprintf('Coincidence Histogram [binwidth = %i ms]',coinsBinWidth))

subplot(2,4,3)
% always -50:50 for jeromiah_jpsth
ii = find(jpsth.xCorrHist(:,1)==-jer_default_lag);
jj = find(jpsth.xCorrHist(:,1)==jer_default_lag);
%bar(jpsth.xCorrHist(ii:jj,1),jpsth.xCorrHist(ii:jj,2))
area(jpsth.xCorrHist(ii:jj,1),smoothdata(jpsth.xCorrHist(ii:jj,2),'gaussian' ,25,'omitnan'),'EdgeColor','none')
title('Cross Correlation Histogram')
subplot(2,4,4)
text(0.25,0.5,'Not Yet')
title('Brody Covariogram - not yet')

subplot(2,4,5)
imagesc(flipud(jer.normalizedJPSTH))
title('Jer-Normalized JPSTH')
subplot(2,4,6)
area(jpsth.coincidenceHist(:,1),smoothdata(jpsth.coincidenceHist(:,2),'gaussian' ,25,'omitnan'),'EdgeColor','none')
title(sprintf('Jer-Coincidence Histogram [binwidth = %i ms]',coinsBinWidth))
subplot(2,4,7)
% always -50:50 for jeromiah_jpsth
area(-jer_default_lag:jer_default_lag,smoothdata(jer.xcorrHist,'gaussian' ,25,'omitnan'),'EdgeColor','none')
title('Jer-Cross Correlation Histogram')

subplot(2,4,8)
area(-jer_default_lag:jer_default_lag,smoothdata(jer.xcorrHist,'gaussian' ,25,'omitnan'),'EdgeColor','none')
hold on
plot(-jer_default_lag:jer_default_lag,jer.sigHigh,'--b')
plot(-jer_default_lag:jer_default_lag,jer.sigLow,'--b')
title(sprintf('Jer-Brody Covariogram [lag = %i ms]',jer_default_lag))
hold off

%%
figure
set(gcf,'Units','normalized')
imagesc(flipud(jpsth.normalizedJpsth),[-1 1]);
set(gca,'YDir','normal')
set(gca,'XTick',[-100:100:500],'YTick',[-100:100:500])
set(gca,'XTickLabel',[],'YTickLabel',[])

colorbar('northoutside')

figure
h_psthx = subplot(2,6,1);
plot(jpsth.xPsthBins{1},jpsth.xPsth{1})
xlim([-100 400])
set(gca,'YDir','reverse')

h_psth6 = subplot(2,6,2);
plot(jpsth.yPsth{1},jpsth.yPsthBins{1})
ylim([-100 400])
set(gca,'XDir','reverse')


