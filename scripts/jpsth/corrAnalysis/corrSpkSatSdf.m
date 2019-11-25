% plot SDFs for different epochs of units used in spike count correlation
% for paired neurons

% Use static windows
spkCorrDir = 'dataProcessed/analysis/11-18-2019/spkCorr';
fn = fullfile(spkCorrDir,'summary/spkCorrAllPairsStaticNew.mat');
spkCorr = load(fn);
outPdfDir = 'dataProcessed/analysis/11-18-2019/spkCorr/summary/sdf';
if ~exist(outPdfDir,'dir')
    mkdir(outPdfDir);
end
outPdfFn = fullfile(outPdfDir,'satSdf_');

%% Fiter data for use too plot sdfs
% Use following area pairs (these are fields of spkCorr loaded above)
pairAreas = {'SEF_SEF','SEF_FEF','SEF_SC'}';
spkCorrBase = cellfun(@(x) fullfile(spkCorrDir,['spkCorr_' x],'mat'),pairAreas,'UniformOutput',false);
% Use following fields
useCols = {
    'srcFile'                      
    'pairAreas'
    'condition'
    'alignedName'
    'alignedEvent'
    'nTrials'                      
    'Pair_UID'                     
    'X_monkey'                     
    'X_sess'                       
    'X_unitNum'                    
    'Y_unitNum'
    'X_area'
    'Y_area'
    'rhoRaw_200ms'                 
    'pvalRaw_200ms'                
 };
%%
spkCorrAllTbl = table();
for pa = 1:numel(pairAreas)
    unitsPlus = spkCorr.(pairAreas{pa})(:,useCols);
    unitsPlus.srcFile = strcat(spkCorrBase{pa},filesep,unitsPlus.srcFile);
%     temp.isSignifByEpoch = temp.pvalRaw_200ms <= useSignif;
%     temp.isSignifPair = cellfun(@(x) sum(temp.isSignifByEpoch(strcmp(temp.Pair_UID,x))) > 0,temp.Pair_UID);
    spkCorrAllTbl = [spkCorrAllTbl; unitsPlus];
    clear temp
end

%% Filter criteria
% Use epoch
useEpoch = 'PostSaccade';
% Use significance level
useSignif = 0.01;
% Use rho percentile
useRhoPercentile = 90;

isEpoch = strcmp(spkCorrAllTbl.alignedName,useEpoch);
isSignifIdx = spkCorrAllTbl.pvalRaw_200ms <= useSignif;
plusIdx = spkCorrAllTbl.rhoRaw_200ms>=0;
minusIdx = spkCorrAllTbl.rhoRaw_200ms<0;

hiPlus = prctile(spkCorrAllTbl.rhoRaw_200ms(plusIdx & isEpoch),useRhoPercentile);
loMinus = prctile(spkCorrAllTbl.rhoRaw_200ms(minusIdx & isEpoch),useRhoPercentile);

hiPlusSignifIdx = spkCorrAllTbl.rhoRaw_200ms>=hiPlus & isSignifIdx & isEpoch;
loMinusSignifIdx = spkCorrAllTbl.rhoRaw_200ms<=loMinus & isSignifIdx & isEpoch;


%% Get units list for + and - spike corrs

plusTable = spkCorrAllTbl(hiPlusSignifIdx,:);
minusTable = spkCorrAllTbl(loMinusSignifIdx,:);
% plus spk corr table

unitsPlus = table();
unitsPlus.unitNum = [plusTable.X_unitNum;plusTable.Y_unitNum];
unitsPlus.area = [plusTable.X_area;plusTable.Y_area];
unitsPlus.threshPlus = repmat(hiPlus,size(unitsPlus,1),1);
unitsPlus = sortrows(unique(unitsPlus,'stable'),'area');
% minus spk corr table
unitsMinus = table();
unitsMinus.unitNum = [minusTable.X_unitNum;minusTable.Y_unitNum];
unitsMinus.area = [minusTable.X_area;minusTable.Y_area];
unitsMinus.threshMinus = repmat(loMinus,size(unitsMinus,1),1);
unitsMinus = sortrows(unique(unitsMinus,'stable'),'area');

%%
unitTbl = table();
uniqUnits = unique([unitsPlus.unitNum;unitsMinus.unitNum]);

for ii = 1:numel(uniqUnits)
    unitNum = uniqUnits(ii);
    pIdx = max([find(unitsPlus.unitNum==unitNum),0]);
    mIdx = max([find(unitsMinus.unitNum==unitNum),0]);
    temp = table();
    if pIdx && ~mIdx
        temp.unitNum = unitNum;
        temp.area = unitsPlus.area(pIdx);
        temp.threshPlus = unitsPlus.threshPlus(pIdx);
        temp.threshMinus = NaN;
        temp.isPlus = 1;
        temp.isMinus = 0;
        temp.isBoth = 0;
    elseif mIdx && ~pIdx
        temp.unitNum = unitNum;
        temp.area = unitsMinus.area(mIdx);
        temp.threshPlus = NaN;
        temp.threshMinus = unitsMinus.threshMinus(mIdx); 
        temp.isPlus = 0;
        temp.isMinus = 1;
        temp.isBoth = 0;
    elseif pIdx && mIdx
        temp.unitNum = unitNum;
        temp.area = unitsMinus.area(pIdx);
        temp.threshPlus = unitsPlus.threshPlus(pIdx);
        temp.threshMinus = unitsMinus.threshMinus(mIdx);  
        temp.isPlus = 0;
        temp.isMinus = 0;
        temp.isBoth = 1;
    end
    unitTbl = [unitTbl;temp];
end

Z = sortrows(unitTbl,{'area','isPlus','isMinus','isBoth','unitNum'});

