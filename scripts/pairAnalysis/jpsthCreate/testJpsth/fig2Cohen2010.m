% Visual : fecfef_m260 1a(RF 3), 2a(RF 3) (target onset: -50 244)
% Movement : fecfef_m250 3b(RF 3,4,5), 4b(RF 5, 6) (target onset : -50 315)

vis = load('data/jpsthFig2/fecfef_m260','-mat');
vxUnitName = 'DSP01a'; vyUnitName = 'DSP02a';
vxCell = vis.(vxUnitName); vyCell = vis.(vyUnitName);
vTargetLoc = vis.P_TargetOn(:,2);
vTargetOn = vis.M_TargetOn;% MNAP target on
vCorrectVar = vis.Correct_(:,2);

allTrials = sum(~isnan(vTargetOn));

vTimeWin = [-50 250];
vCorrectTrialsRFAll = find(vCorrectVar==1);
vCorrectTrialsRF = find(vCorrectVar==1 & vTargetLoc==3);

selTrials = vCorrectTrialsRFAll;

vXAligned = SpikeUtils.alignSpikeTimes(vxCell(selTrials,:),vTargetOn(selTrials),vTimeWin);
vYAligned = SpikeUtils.alignSpikeTimes(vyCell(selTrials,:),vTargetOn(selTrials),vTimeWin);
binWidth = 1;
jpsth = SpikeUtils.jpsth(vXAligned,vYAligned,vTimeWin,binWidth,50);

jer = SpikeUtils.jeromiahJpsth(vXAligned,vYAligned,vTimeWin,binWidth,50);


xCounts = [ max(sum(jpsth.xSpikeCounts{1})) max(sum(jer.xCounts))]
yCounts = [ max(sum(jpsth.ySpikeCounts{1})) max(sum(jer.yCounts))]

xPsthMax = ([max(jer.psth_1) max(jpsth.xPsth{1})]).*1000/1
yPsthMax = ([max(jer.psth_2) max(jpsth.yPsth{1})]).*1000/1


