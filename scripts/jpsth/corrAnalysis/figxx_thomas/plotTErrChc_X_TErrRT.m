function [ ] = plotTErrChc_X_TErrRT( ninfo , nstats )
%plotTErrChc_X_TErrRT Summary of this function goes here
%   Detailed explanation goes here

idxArea = ismember({ninfo.area}, {'SEF'});
idxDa = ismember({ninfo.monkey}, {'D'});
idxEu = ismember({ninfo.monkey}, {'E'});

idxErr = ([ninfo.errGrade] >= 2);
idxRew = (abs([ninfo.rewGrade]) >= 2 & ~isnan([nstats.A_Reward_tErrStart_Fast]));

idxBothErrDa = (idxArea & idxDa & (idxErr & idxRew));
idxBothErrEu = (idxArea & idxEu & (idxErr & idxRew));

%compute median time of error signaling for each monkey
tmedChcErrDa = median([nstats(idxErr & idxDa).A_ChcErr_tErr_Fast]);
tmedChcErrEu = median([nstats(idxErr & idxEu).A_ChcErr_tErr_Fast]);
tmedRTErrDa = median([nstats(idxRew & idxDa).A_Reward_tErrStart_Acc]);
tmedRTErrEu = median([nstats(idxRew & idxEu).A_Reward_tErrStart_Acc]);

%collect latencies for neurons with both types of modulation
tChcErrDa = [nstats(idxBothErrDa).A_ChcErr_tErr_Fast];% - tmedChcErrDa;
tChcErrEu = [nstats(idxBothErrEu).A_ChcErr_tErr_Fast];% - tmedChcErrEu;
tRTErrDa = [nstats(idxBothErrDa).A_Reward_tErrStart_Acc];% - tmedRTErrDa;
tRTErrEu = [nstats(idxBothErrEu).A_Reward_tErrStart_Acc];% - tmedRTErrEu;

%combine corrected latencies across monkeys
tChcErr = [tChcErrDa , tChcErrEu];
tRTErr  = [tRTErrDa , tRTErrEu];


%% Plotting

figure(); hold on
scatter(tChcErr, tRTErr, 30, 'k', 'filled')
ppretty([4.8,3])

end%fxn:plotTErrChc_X_TErrRT()

