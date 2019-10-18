function [ CRtbl] = plotSDF_ErrorCompare_SAT( binfo , moves , ninfo , nstats , spikes , varargin )
%plotSDF_ErrorCompare_SAT() Summary of this function goes here
%   Detailed explanation goes here

%args = getopt(varargin, {{'area=','SEF'}, {'monkey=',{'D','E','Q','S'}}});
args = getopt(varargin, {{'area=','SEF'}, {'monkey=',{'D','E'}}});

idxArea = ismember({ninfo.area}, args.area);
idxMonkey = ismember({ninfo.monkey}, args.monkey);

idxErr = ([ninfo.errGrade] >= 2);
idxRew = (abs([ninfo.rewGrade]) >= 2 & ~isnan([nstats.A_Reward_tErrStart_Fast]));

idxKeep = (idxArea & idxMonkey & (idxErr | idxRew));

NUM_CELLS = sum(idxKeep);
ninfo = ninfo(idxKeep);
spikes = spikes(idxKeep);

tplotResp = (-200 : 450);   offsetResp = 200;
tplotRew  = (-50 : 600);    offsetRew = 50;

%initialization -- contrast ratio (A_err - A_corr) / (A_err + A_corr)
CRacc = NaN(1,NUM_CELLS);   tEstAcc = (100 : 500) + offsetRew;
CRfast = NaN(1,NUM_CELLS);  tEstFast = (100 : 300) + offsetResp;
CRunit = NaN(1,NUM_CELLS);
CRmonk = cell(1,NUM_CELLS);
for cc = 1:NUM_CELLS
  fprintf('%s - %s - %d\n', ninfo(cc).sess, ninfo(cc).unit, ninfo(cc).unitNum)
  CRunit(cc) = ninfo(cc).unitNum;
  CRmonk{cc} = ninfo(cc).monkey;
  kk = ismember({binfo.session}, ninfo(cc).sess);
  
  rtKK = double(moves(kk).resptime);
  trewKK = double(binfo(kk).rewtime + binfo(kk).resptime);
  
  %index by isolation quality
  idxIso = identify_trials_poor_isolation_SAT(ninfo(cc), binfo(kk).num_trials);
  %index by condition
  idxAcc = (binfo(kk).condition == 1 & ~idxIso & ~isnan(trewKK));
  idxFast = (binfo(kk).condition == 3 & ~idxIso & ~isnan(trewKK));
  %index by trial outcome
  idxCorr = ~(binfo(kk).err_dir | binfo(kk).err_time | binfo(kk).err_hold | binfo(kk).err_nosacc);
  idxErrChc = (binfo(kk).err_dir & ~binfo(kk).err_time);
  idxErrTime = (~binfo(kk).err_dir & binfo(kk).err_time);
  %index by screen clear on Fast trials
  idxClear = logical(binfo(kk).clearDisplayFast);
  
  %compute single-trial SDF
  sdfStim = compute_spike_density_fxn(spikes(cc).SAT);
  sdfResp = align_signal_on_response(sdfStim, rtKK);
  sdfRew  = align_signal_on_response(sdfStim, trewKK);
  
  %split SDF into groups and compute mean
  sdfFastCorr = nanmean(sdfResp(idxFast & idxCorr, tplotResp + 3500));
  sdfFastErr = nanmean(sdfResp(idxFast & idxErrChc & ~idxClear, tplotResp + 3500));
  sdfAccCorr = nanmean(sdfRew(idxAcc & idxCorr, tplotRew + 3500));
  sdfAccErr  = nanmean(sdfRew(idxAcc & idxErrTime, tplotRew + 3500));
  
  %compute contrast ratio
  muAccCorr = mean(sdfAccCorr(tEstAcc));      muAccErr = mean(sdfAccErr(tEstAcc));
  muFastCorr = mean(sdfFastCorr(tEstFast));   muFastErr = mean(sdfFastErr(tEstFast));
  
  CRacc(cc) = (muAccErr - muAccCorr) / (muAccErr + muAccCorr);
  CRfast(cc) = (muFastErr - muFastCorr) / (muFastErr + muFastCorr);
  
end%for:cells(cc)

%% Plotting - Contrast ratio

figure(); hold on
scatter(CRfast, CRacc, 30, 'r', 'filled')
xlabel('Contrast ratio - Choice error')
ylabel('Contrast ratio - Timing error')
%ppretty([4.8,3])

CRtbl = table();
CRtbl.monk = CRmonk';
CRtbl.unitNum = CRunit';
CRtbl.CRfast = CRacc';
CRtbl.CRacc = CRacc';


end%fxn:plotSDF_ErrorCompare_SAT()
