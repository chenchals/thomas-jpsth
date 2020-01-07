function [ CRtbl] = plotSDF_ErrorCompare_SAT( binfo , moves , unitInfo , unitStats , spikes , varargin )
%plotSDF_ErrorCompare_SAT() Summary of this function goes here
%   Detailed explanation goes here

%args = getopt(varargin, {{'area=','SEF'}, {'monkey=',{'D','E','Q','S'}}});
args = getopt(varargin, {{'area=','SEF'}, {'monkey=',{'D','E'}}});

useRewGradeFilter = 0;
% in place of [nstats.A_Reward_tErrStart_Fast] old
% use unitStats.ErrorChoiceSignalTimFastStart
errorChoiceSignalTimFastStart = unitStats.ErrorChoiceSignalTimFastStart;

idxArea = ismember(unitInfo.area, args.area);
idxMonkey = ismember(unitInfo.monkey, args.monkey);

idxAllErr = abs(unitInfo.errGrade) >= 2;
idxAllRew = abs(unitInfo.rewGrade) >= 2 ; %& ~isnan(errorChoiceSignalTimFastStart);

if useRewGradeFilter
    idxKeep = (idxArea & idxMonkey & (idxAllErr | idxAllRew) & ~isnan(errorChoiceSignalTimFastStart));
else
    idxKeep = (idxArea & idxMonkey & (idxAllErr | idxAllRew));
end

idxKeep = (idxArea & idxMonkey & (idxAllErr | idxAllRew));

NUM_CELLS = sum(idxKeep);
unitInfo = unitInfo(idxKeep,:);
spikes = spikes(idxKeep);

tplotResp = (-200 : 450);   offsetResp = 200;
tplotRew  = (-50 : 600);    offsetRew = 50;

%initialization -- contrast ratio (A_err - A_corr) / (A_err + A_corr)
CRacc = nan(NUM_CELLS,1);   tEstAcc = (100 : 500) + offsetRew;
CRfast = nan(NUM_CELLS,1);  tEstFast = (100 : 300) + offsetResp;
CRunit = nan(NUM_CELLS,1);
CRmonk = cell(NUM_CELLS,1);
CRsession = cell(NUM_CELLS,1);
CRerrGrade = nan(NUM_CELLS,1);
CRrewGrade = nan(NUM_CELLS,1);
% another filter for rewGrade: ~isnan([nstats.A_Reward_tErrStart_Fast])
CRrewGrade_A_Reward_tErrStart_Fast = nan(NUM_CELLS,1);
CRcount = nan(NUM_CELLS,1);
CRisErrOnly = zeros(NUM_CELLS,1);
CRisRewOnly = zeros(NUM_CELLS,1);
CRisBothErrRew = zeros(NUM_CELLS,1);

for cc = 1:NUM_CELLS
  fprintf('#%d : %s - %s - %d\n', cc, unitInfo.sess{cc}, unitInfo.unit{cc}, unitInfo.unitNum(cc))
  kk = ismember(binfo.session, unitInfo.sess{cc});
  
  rtKK = double(moves.resptime{kk});
  trewKK = double(binfo.rewtime{kk} + binfo.resptime{kk});
  
  %index by isolation quality
  idxIso = identify_trials_poor_isolation_SAT(unitInfo(cc,:), binfo.num_trials(kk));
  %index by condition
  idxAcc = (binfo.condition{kk} == 1 & ~idxIso & ~isnan(trewKK));
  idxFast = (binfo.condition{kk} == 3 & ~idxIso & ~isnan(trewKK));
  %index by trial outcome
  idxCorr = ~(binfo.err_dir{kk} | binfo.err_time{kk} | binfo.err_hold{kk} | binfo.err_nosacc{kk});
  idxErrChc = (binfo.err_dir{kk} & ~binfo.err_time{kk});
  idxErrTime = (~binfo.err_dir{kk} & binfo.err_time{kk});
  %index by screen clear on Fast trials
  idxClear = logical(binfo.clearDisplayFast{kk});
  
  %compute single-trial SDF
  sdfStim = compute_spike_density_fxn(spikes{cc});
  sdfResp = align_signal_on_response(sdfStim, rtKK);
  sdfRew  = align_signal_on_response(sdfStim, trewKK);
  
  %split SDF into groups and compute mean
  sdfFastCorr = nanmean(sdfResp(idxFast & idxCorr, tplotResp + 3500));
  sdfFastErr = nanmean(sdfResp(idxFast & idxErrChc & ~idxClear, tplotResp + 3500));
  sdfAccCorr = nanmean(sdfRew(idxAcc & idxCorr, tplotRew + 3500));
  sdfAccErr  = nanmean(sdfRew(idxAcc & idxErrTime, tplotRew + 3500));
  
  %compute contrast ratio
  muAccCorr = mean(sdfAccCorr(tEstAcc));      
  muAccErr = mean(sdfAccErr(tEstAcc));
  muFastCorr = mean(sdfFastCorr(tEstFast));   
  muFastErr = mean(sdfFastErr(tEstFast));
  
  CRacc(cc) = (muAccErr - muAccCorr) / (muAccErr + muAccCorr);
  CRfast(cc) = (muFastErr - muFastCorr) / (muFastErr + muFastCorr);
  
  % additional criteria for plot 
  isErr = abs(unitInfo.errGrade(cc)) >=2 & ~(abs(unitInfo.rewGrade(cc)) >=2);
  isRew = ~(abs(unitInfo.errGrade(cc)) >=2) & abs(unitInfo.rewGrade(cc)) >=2;
  isBoth = abs(unitInfo.errGrade(cc)) >=2 & abs(unitInfo.rewGrade(cc)) >=2;
  CRisErrOnly(cc) = isErr;
  CRisRewOnly(cc) = isRew;
  CRisBothErrRew(cc) = isBoth;
    
  CRcount(cc) = cc;
  CRunit(cc) = unitInfo.unitNum(cc);
  CRmonk{cc} = unitInfo.monkey{cc};
  CRsession{cc} = unitInfo.sess{cc};
  CRerrGrade(cc) = unitInfo.errGrade(cc);
  CRrewGrade(cc) = unitInfo.rewGrade(cc);
  % another filter for rewGrade: ~isnan([nstats.A_Reward_tErrStart_Fast])
  CRrewGrade_A_Reward_tErrStart_Fast(cc) = errorChoiceSignalTimFastStart(cc);
  
  
end%for:cells(cc)
%% out table
CRtbl = table();
CRtbl.counter = CRcount;
CRtbl.monk = CRmonk;
CRtbl.session = CRsession;
CRtbl.unitNum = CRunit;
CRtbl.CRerrGrade =  CRerrGrade;
CRtbl.CRrewGrade =  CRrewGrade ;
% another filter for rewGrade: ~isnan([nstats.A_Reward_tErrStart_Fast])
CRtbl.filter_rewGrade_A_Reward_tErrStart_Fast =  CRrewGrade_A_Reward_tErrStart_Fast;
CRtbl.CRfast = CRfast;
CRtbl.CRacc = CRacc;
CRtbl.CRisErrOnly = CRisErrOnly;
CRtbl.CRisRewOnly = CRisRewOnly;
CRtbl.CRisBothErrRew = CRisBothErrRew;
% Add quadrant num to table
qx = CRtbl.CRfast; qy = CRtbl.CRacc;
CRtbl.quadrantNum(qx >= 0 & qy >= 0) = 1;
CRtbl.quadrantNum(qx < 0 & qy >= 0) = 2;
CRtbl.quadrantNum(qx < 0 & qy < 0) = 3;
CRtbl.quadrantNum(qx >= 0 & qy < 0) = 4;


%% Plotting - Contrast ratio

figure(); hold on
legTxt = {};
%scatter(CRfast,CRacc,50,'s');
% % is Err only
ce = CRtbl.CRisErrOnly==1;
g = grpstats(CRtbl(ce,{'monk','CRisErrOnly'}),'monk');
t = arrayfun(@(x) sprintf('%s:%d',char(g{x,1}),g{x,2}),1:2,'UniformOutput',false);
legTxt = [legTxt,['Error only (' char(join(t,', ')) ')']];
scatter(CRtbl.CRfast(ce),CRtbl.CRacc(ce),80,'filled','r','MarkerFaceColor',[0.5 0.2 0.0]);
% % is Rew only
te = CRtbl.CRisRewOnly==1;
g = grpstats(CRtbl(te,{'monk','CRisRewOnly'}),'monk');
t = arrayfun(@(x) sprintf('%s:%d',char(g{x,1}),g{x,2}),1:2,'UniformOutput',false);
legTxt = [legTxt,['Reward only (' char(join(t,', ')) ')']];
scatter(CRtbl.CRfast(te),CRtbl.CRacc(te),80,'filled','^','MarkerFaceColor',[0.0 0.5 0.2]);
% % is both Err and Rew 
cte = CRtbl.CRisBothErrRew==1;
g = grpstats(CRtbl(cte,{'monk','CRisBothErrRew'}),'monk');
t = arrayfun(@(x) sprintf('%s:%d',char(g{x,1}),g{x,2}),1:2,'UniformOutput',false);
legTxt = [legTxt,['Error & Reward (' char(join(t,', ')) ')']];
scatter(CRtbl.CRfast(cte),CRtbl.CRacc(cte),80,'filled','s','MarkerFaceColor',[0.2 0.2 0.5]);

xlim(round(minmax(CRtbl.CRfast')+[-0.05 0.05],2))
ylim(round(minmax(CRtbl.CRacc')+[-0.05 0.05],2))

xlabel('Contrast ratio - Choice error')
ylabel('Contrast ratio - Timing error')

line([0 0],get(gca,'YLim'),'color','k','HandleVisibility','off');
line(get(gca,'XLim'),[0 0],'color','k','HandleVisibility','off');

legend(legTxt,'Location','northeast','Box','off','Orientation','vertical')

ppretty([4.8,3])


end%fxn:plotSDF_ErrorCompare_SAT()
