% Compute and plot SDFs for every outcome split on outcome of previous trial
% cOC: current Outcome
% pOC: previous outcome [1,2 or 3]
% A : Accurate,  F: Fast
% _C : Correct, _EC: ErrorChoice, _ET: ErrorTiming
% _____________________________
% | cOC  | pOC1 | pOC2 | pOC3 |  
% |------|------|------|------|
% | A_C  | A_C  | A_EC | A_ET |
% | A_EC | A_C  | A_EC | A_ET |
% | A_ET | A_C  | A_EC | A_ET |
% |------|------|------|------|
% | F_C  | F_C  | F_EC | F_ET |
% | F_EC | F_C  | F_EC | F_ET |
% | F_ET | F_C  | F_EC | F_ET |
% -----------------------------
%

trialHistoryFile = 'dataProcessed/analysis/spkCorr/trialHistory/SpkCorr_TrialHistory.mat';
trialHist = load(trialHistoryFile);
% date for trialHistory for previousTrial = trialHistory_1
trialHist = trialHist.trialHistory_1;
% Rows are: 
% nOutcomes (6: 3Accu, 3Fast) * nSAT-Sessions-With-Pairs (14: 8D, 6E)
% Columns are:
% session, currOutcome, currOutcomeTrls, currOutcomeTrlsCount, 
%   *_p1_count: count of number of trials for previous (-1) outcome
% FastCorrect_p1_count, FastErrorChoice_p1_count, FastErrorTiming_p1_count,
% AccurateCorrect_p1_count, AccurateErrorChoice_p1_count,
% AccurateErrorTiming_p1_count, 
%   *_p1: trial nos for previous (-1) outcome
%   *_p1 + 1 = current trial outcome (currOutcome column) 
%              given previous trial outcome
% FastCorrect_p1, FastErrorChoice_p1, FastErrorTiming_p1,
% AccurateCorrect_p1, AccurateErrorChoice_p1, AccurateErrorTiming_p1  
colNames = trialHist.Properties.VariableNames';

%% Plot num. trials previous outcome for each outcome group unaccounted curr Trials under other-outcome
% That is, nCurrOutcome ~= prevCorr + prevErrChoice + prevErrTiming
currOutcomes = {'Correct','ErrorChoice','ErrorTiming'}';
prevOutcomes = [currOutcomes;'Other'];
fastColNamesCount = colNames(~cellfun(@isempty,regexp(colNames,'Fast.*[Cc]ount$','match')));
accuColNamesCount = colNames(~cellfun(@isempty,regexp(colNames,'Accurate.*[Cc]ount$','match')));

fastCols = ['currOutcome';'currOutcomeTrlsCount';fastColNamesCount];
accuCols = ['currOutcome';'currOutcomeTrlsCount';accuColNamesCount];

fastTbl = grpstats(trialHist(contains(trialHist.currOutcome,'Fast'),fastCols),'currOutcome','sum');
fastTbl.sum_FastOther_p1_count = fastTbl{:,3} - sum(fastTbl{:,4:6},2);
fastData = array2table(fastTbl{:,4:7}./fastTbl{:,3},'VariableNames',prevOutcomes);
fastData.outcome = currOutcomes;

accuTbl = grpstats(trialHist(contains(trialHist.currOutcome,'Accu'),accuCols),'currOutcome','sum');
accuTbl.sum_AccurateOther_p1_count = accuTbl{:,3} - sum(accuTbl{:,4:6},2);

accuData = array2table(accuTbl{:,4:7}./accuTbl{:,3},'VariableNames',prevOutcomes);
accuData.outcome = currOutcomes;

H_fig = newFigure();
accColor = [1 0 0];
fastColor = [0 0.5 0];
temp = [fastData{:,prevOutcomes};accuData{:,prevOutcomes}];
yLim = [0 round(max(temp(:))+0.005,2)];
outcomes = prevOutcomes; % includes 'other'
for oc = 1:numel(currOutcomes)
    outcome = outcomes{oc};
    yLab = sprintf('Prob.[prevTrialOutcome|\ncurrTrialOutcome = %s]',outcome);
    titleStr = sprintf('Current trial outcome : %s',outcome);
    subplot(3,1,oc);
    yDat = [fastData{oc,prevOutcomes};accuData{oc,prevOutcomes}]';
    h = bar(1:numel(prevOutcomes), yDat);
    set(h(1),'FaceColor',fastColor);set(h(2),'FaceColor',accColor);
    set(h,'FaceAlpha',0.5)
    set(gca,'xticklabels',prevOutcomes,'fontWeight','bold')  
    ylabel(yLab)
    ylim(yLim);
    title(titleStr); 
    legend({'Fast','Accurate'},'Location','north','Orientation','horizontal','Box','off')
end

saveFigPdf('test.pdf')




