% compute and compare JPSTH
fDir = '/Volumes/schalllab/data/Fechner/fecfef/MatlabDataFiles';

sess = 'fecfef_m226.mat';

fdata = load(fullfile(fDir,sess),'-mat');
xUnitName = 'DSP09a'; yUnitName = 'DSP09c';
xCell = fdata.(xUnitName); yCell = fdata.(yUnitName);
targetLoc = fdata.P_TargetOn(:,2);
targetOn = fdata.M_TargetOn;% MNAP target on
correctVar = fdata.Correct_(:,2);

allTrials = sum(~isnan(targetOn));

timeWin = [0 300];

correctTrialsRFAll = find(correctVar==1);
correctTrialsRF = find(correctVar==1 & (targetLoc==3 | targetLoc==7));

selTrials = correctTrialsRF;

xAligned = SpikeUtils.alignSpikeTimes(xCell(selTrials,:),targetOn(selTrials),timeWin);
yAligned = SpikeUtils.alignSpikeTimes(yCell(selTrials,:),targetOn(selTrials),timeWin);
binWidth = 5;
jpsth = SpikeUtils.jpsth(xAligned,yAligned,timeWin,binWidth,5);

jer = SpikeUtils.jeromiahJpsth(xAligned,yAligned,timeWin,binWidth,5);


xCounts = [ max(sum(jpsth.xSpikeCounts{1})) max(sum(jer.xCounts))]
yCounts = [ max(sum(jpsth.ySpikeCounts{1})) max(sum(jer.yCounts))]

xPsthMax = ([max(jer.psth_1) max(jpsth.xPsth{1})]).*1000/1
yPsthMax = ([max(jer.psth_2) max(jpsth.yPsth{1})]).*1000/1

