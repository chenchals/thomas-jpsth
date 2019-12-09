function [spkCorrSdfs] = corrSpkSatSdf(outcomeToUse,epochToUse)
% CORRSPKSATSDF Compute SDFs using filter criteria
%    outcomeToUse : valid values are :'Correct', 'ErrorChoice',
%                   'ErrorTiming' 
%    epochToUse : valid values are :'Baseline', 'Visual', 'PostSaccade',
%                 'PostReward' 
%  Useage: 
%     corrSpkSatSdf('Correct','Baseline'): for Correct trials and Baseline
%                   epoch, that computes data and saves file
%                   spkCorrSdfs_Correct_Baseline.mat  
%
% Computes SDFs units with significant spike correlations with any other
% unit(s) in the session in the PostSaccade epoch. Use data from:
% analysis/11-18-2019/spkCorr/summary/spkCorrAllPairsStaticNew.mat
% dataset/spikes_SAT.mat
% dataset/dataNeurophys_SAT.mat
% dataset/TrialTypesDB.mat
% dataset/TrialEventTimesDB.mat
% see also CORRSPKSTATICEXTRACT

spkCorrDir = 'dataProcessed/analysis/11-18-2019/spkCorr';
% for loading event times and trial types
datasetDir = 'dataProcessed/dataset';
% Files from which to load spik corrs and compute SDFs
spkCorrFile = fullfile(spkCorrDir,'summary/spkCorrAllPairsStaticNew.mat');
spikeTimesFile = fullfile(datasetDir,'spikes_SAT.mat');
unitInfoStatsFile = fullfile(datasetDir,'dataNeurophys_SAT.mat');
trialTypesFile = fullfile(datasetDir,'TrialTypesDB.mat');
trialEventTimesFile = fullfile(datasetDir,'TrialEventTimesDB.mat');

%% Filter criteria - to limit units that will be processed for SDFs
% available conditions

%conditions
availableConditions = {
    'AccurateCorrect';'AccurateErrorChoice';'AccurateErrorTiming';
    'FastCorrect';    'FastErrorChoice';    'FastErrorTiming'
    };
% Use epoch available: {Baseline, Visual, PostDSaccade, PostReward}
% These are column names / values in the spkCorrAllPairsStaticNew.mat
useEpoch = epochToUse; %'PostSaccade'; % value for alignedName column
% Use condition that contain the outcome values below
useOutcome = outcomeToUse;%'ErrorTiming'; % condition column values shall contain this string
% Output specs - add epoch and outcome used
spkCorrSdfsFile = fullfile(spkCorrDir,['summary/spkCorrSdfs_' useOutcome '_' useEpoch '.mat']);
% Use unsigned rho
useRhoUnsigned = true;
% Use rho value from column
useRhoValColName = 'rhoRaw_150ms';
% Use p value from column
usePvalColName = 'pvalRaw_150ms';
% Use rho percentile
useRhoPercentile = 70;
% Use significance level for selecting units to plot
useSignif = 0.01;
% For a given unit, use significance level for selecting paired units 
usePairSignif = 0.05;
% Selected trials for any condition must be greater than this threshold
useMinTrialCount = 50;
% use Conditions
% useConditions = availableConditions(contains(availableConditions,useOutcome));
useConditions = availableConditions; % so we can compute SDFs for all conditions


%% Parameters for SDFs
useAreas = {'SEF' 'FEF' 'SC'};
% Setup time windows for different event time alignment, the field names
% SHALL correspond to column names for trialEventTimes below.
alignEvents = {'CueOn','SaccadePrimary','RewardTime'};
alignTimeWins = {[-600 400],[-200 600],[-100 400]};
alignNames = {'Visual','PostSaccade','PostReward'};
% first sort by accurate/fast?
firstSortEventNames = {'SaccadePrimary','SaccadePrimary','SaccadePrimary'};
% use epsp kernel length to pad alignTimeWin
% pad one-half the length of epsp kernel
% In convn call, use flag 'valid' to drop off the padded window
pspKernel = getPspKernelForSdf();
% pspKernel must be even
if ~mod(numel(pspKernel),2)
    pspKernel(end+1) = pspKernel(end);% last point
end
% Pad alignTimeWin
padLen = floor(numel(pspKernel)/2);
padTimeWin = [-padLen padLen];

% Use following area pairs 
% These are variable names in spkCorrAllPairsStaticNew.mat
usePairAreas = {'SEF_SEF','SEF_FEF','SEF_SC'}';
% Use following column names from each pairArea table
% These are column names in variables in spkCorrAllPairsStaticNew.mat
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
    %'Y_sess' for a pair session is always the same
    'X_unitNum'
    'Y_unitNum'
    'X_area'
    'Y_area'
    useRhoValColName
    usePvalColName
    };

%% Load data variables
stTime = tic;
fprintf('Loading data variables...');

if ~exist('spkCorr','var')
    % Load data variable: spike correlations
    spkCorr = load(spkCorrFile);
    % Load data variable: spike times
    spikesSat = load(spikeTimesFile);
    spikesSat = spikesSat.spikesSAT;
    % Load unit Info for knowing which trials to be removed of any
    unitInfoAll = load(unitInfoStatsFile,'unitInfo');
    unitInfoAll = unitInfoAll.unitInfo;
    % Load data variable: TrialTypes
    sessionTrialTypes = load(trialTypesFile);
    sessionTrialTypes = sessionTrialTypes.TrialTypesDB;
    % Load data variable: TrialEventTimes
    sessionEventTimes = load(trialEventTimesFile);
    sessionEventTimes = sessionEventTimes.TrialEventTimesDB;
end
lap = toc(stTime);
fprintf('done %5.2f sec\n',lap);

%% Aggregate data into a single table for areas of intertest
spkCorrAll = table();
for pa = 1:numel(usePairAreas)
    temp =  spkCorr.(usePairAreas{pa})(:,useCols);
    % recode rho, pval
    temp.rho = temp.(useRhoValColName);
    temp.pval = temp.(usePvalColName);
    temp.rhoUnsigned = abs(temp.rho);
    % recode outcome, satCondition
    temp.outcome = regexprep(temp.condition,'(Fast)|(Accurate)','');
    temp.satCondition = regexprep(temp.condition,'(Correct)|(Error.*)','');
    spkCorrAll = [spkCorrAll;temp];
    clearvars temp
end

%% Filter Data to limit Units used for SDFs
isOutcome = strcmp(spkCorrAll.outcome,useOutcome);
isEpoch = strcmp(spkCorrAll.alignedName,useEpoch);
isMinTrialCount = spkCorrAll.nTrials >= useMinTrialCount;
spkCorrAllTbl = spkCorrAll(isOutcome & isEpoch & isMinTrialCount,:);

isSignif = spkCorrAllTbl.pval <= useSignif;


hiUnsigned = prctile(spkCorrAllTbl.rhoUnsigned,useRhoPercentile);
isHiUnsigned = spkCorrAllTbl.rhoUnsigned >= hiUnsigned;
isHiUnsignedSignif = isHiUnsigned & isSignif;

isPlusIdx = spkCorrAllTbl.rho>=0;
hiPlus = prctile(spkCorrAllTbl.rho(isPlusIdx),useRhoPercentile);
isHiPlus = spkCorrAllTbl.rho>=hiPlus;
isHiPlusSignif = spkCorrAllTbl.rho>=hiPlus & isSignif;

isMinusIdx = spkCorrAllTbl.rho<0;
loMinus = prctile(spkCorrAllTbl.rho(isMinusIdx),useRhoPercentile);
isLoMinus = spkCorrAllTbl.rho<=loMinus;
isLoMinusSignif = spkCorrAllTbl.rho<=loMinus & isSignif;


% Filter criteria
filterCriteria = struct();
filterCriteria.usePairAreas = usePairAreas;
filterCriteria.useEpoch = useEpoch;
filterCriteria.useOutcome = useOutcome;
filterCriteria.useRhoPercentile = useRhoPercentile;
filterCriteria.usePvalForUnit = useSignif;
filterCriteria.useMinTrialCount = useMinTrialCount;
filterCriteria.useRhoUnsigned = useRhoUnsigned;

filterCriteria.rhoUnsignedThresh = hiUnsigned;
filterCriteria.rhoPositiveThresh = hiPlus;
filterCriteria.rhoNegativeThresh = loMinus;
filterCriteria.usePvalForPairedUnits = usePairSignif;


%% Get list of units for positive and negative Spike correlations
if useRhoUnsigned
    tempTable = spkCorrAllTbl(isHiUnsignedSignif,:);
    unitTbl = getUniqueUnitTbl(tempTable);
else
    tempTable = spkCorrAllTbl(isHiPlusSignif,:);
    tempTable = [tempTable; spkCorrAllTbl(isLoMinusSignif,:)];
    unitTbl = getUniqueUnitTbl(tempTable);
end

%% Process units for SDF by area
fprintf('Processing units for SDF by area\n');
unitSdfs = {};
opts = table();
roNum = 0;
warning('off');
for jj = 1:numel(useAreas)
    currArea = useAreas{jj};
    currUnitsTbl = unitTbl(strcmp(unitTbl.area,currArea),:);
    nUnits = size(currUnitsTbl,1);
    tic
    fprintf('   Area: %s...',currArea);
    s = '';
    for uu = 1:nUnits
        %%
        
        currUnitTbl = currUnitsTbl(uu,:);
        unitNum = currUnitTbl.unitNum;
        currUnitInfo = unitInfoAll(unitInfoAll.unitNum==unitNum,:);
        
        sess = currUnitTbl.sess{1};
        unit = currUnitInfo.unit{1};
        rho = currUnitTbl.rho{1};
        rhoUnsigned = currUnitTbl.rhoUnsigned{1};
        rscIsPositive = currUnitTbl.isPlus;
        rscIsNegative = currUnitTbl.isMinus;
        rscIsBoth = currUnitTbl.isBoth;      
        
        % fun stuff...
        fprintf(repmat('\b',1,size(s,2)));
        s = sprintf('Unit %d [%s-%s] %d of %d',unitNum,sess,unit,uu,nUnits);
        fprintf('%s',s);
        
        %
        unitSpkTimes = spikesSat{unitNum}';
        evntTimes = sessionEventTimes(strcmp(sessionEventTimes.session,sess),:);
        trialTypes = sessionTrialTypes(strcmp(sessionTrialTypes.session,sess),:);
        trialsToRemove = unitInfoAll.trRemSAT{unitInfoAll.unitNum==unitNum};
        
        % For this unit find spikecorrs that are significant
        isXunit = sum(spkCorrAllTbl.X_unitNum==unitNum)>0;
        if isXunit
            useUnitField = 'X_unitNum';
            usePairUnitField = 'Y_unitNum';
            usePairAreaField = 'Y_area';
        else
            useUnitField = 'Y_unitNum';
            usePairUnitField = 'X_unitNum';
            usePairAreaField = 'X_area';
        end
        signifRscIdx = find(spkCorrAllTbl.(useUnitField)==unitNum & spkCorrAllTbl.pval <= usePairSignif);
        signifUnitTbl = table();
        signifUnitTbl.pairUnitNum = spkCorrAllTbl.(usePairUnitField)(signifRscIdx);
        signifUnitTbl.pairUnitArea = spkCorrAllTbl.(usePairAreaField)(signifRscIdx);
        signifUnitTbl.condition = spkCorrAllTbl.condition(signifRscIdx);
        
        %%
        for cc = 1:numel(useConditions)
            % SDF for each condition,
            % some condition may not exist...so use try...catch
            try
                condition = useConditions{cc};
                % Paired Units that have signif spike corr. for condition
                tempTbl = signifUnitTbl(strcmp(signifUnitTbl.condition,condition),:);
                tempTbl = unique(tempTbl,'stable');
                pairedSefUnitNums = tempTbl.pairUnitNum(strcmp(tempTbl.pairUnitArea,'SEF'));
                pairedFefUnitNums = tempTbl.pairUnitNum(strcmp(tempTbl.pairUnitArea,'FEF'));
                pairedScUnitNums = tempTbl.pairUnitNum(strcmp(tempTbl.pairUnitArea,'SC'));
                
                % Get trials for this unit discounting the trials to remove 
                selTrials = getSelectedTrials(condition,trialTypes,trialsToRemove,useMinTrialCount);
                if isempty(selTrials)
                    continue;
                end
                % Increment row counter for output - pivot alignEvents
                % There is 1 row per condition. The column names capture
                % ALIGNED_NAME_[raters,sdf etc]
                roNum = roNum + 1;
                opts.session{roNum} = sess;
                opts.unit{roNum} = unit;
                opts.unitNum(roNum) = unitNum;
                opts.area{roNum} = currArea;
                opts.condition{roNum} = condition;
                opts.rho{roNum} = rho;
                opts.rhoUnsigned{roNum} = rhoUnsigned;
                opts.rscIsPositive(roNum) = rscIsPositive;
                opts.rscIsNegative(roNum) = rscIsNegative;
                opts.rscIsBoth(roNum) = rscIsBoth;
                % pivot all aligned events. That is make separate columns
                % for each aligned event by prefixing the columnName with
                % aligned name
                for evId = 1:numel(alignEvents)
                    alignedEvent = alignEvents{evId};
                    alignedTimeWin = alignTimeWins{evId}; %#ok<*PFBNS>
                    alignedTimeWinPad = alignedTimeWin + padTimeWin;
                    alignedName = alignNames{evId};
                    % sort trials
                    firstSortByName = firstSortEventNames{evId};
                    if isempty(firstSortByName)
                        firstSortByTime = -inf(numel(selTrials),1);
                    else
                        % sort selected trials based on the sort event name
                        % many ways, but using table to do this
                        temp = double(evntTimes.(firstSortByName){1}(selTrials));
                        temp(temp==0 | isnan(temp)) = -inf;
                        selTrlsMat = [selTrials temp];
                        selTrlsMat = sortrows(selTrlsMat,2);
                        selTrials = selTrlsMat(:,1);
                        firstSortByTime = selTrlsMat(:,2);
                    end
                    % Get align times for sorted trials
                    alignTime = evntTimes.CueOn{1};
                    if ~strcmp(alignedEvent,'CueOn')
                        alignTime = alignTime + double(evntTimes.(alignedEvent){1}(:));
                    end
                    alignTime = alignTime(selTrials);
                    % Align Spike times and get rasters
                    alignedSpkTimes = SpikeUtils.alignSpikeTimes(unitSpkTimes(selTrials),alignTime, alignedTimeWinPad);
                    tempRast = SpikeUtils.rasters(alignedSpkTimes,alignedTimeWinPad);
                    rasters = tempRast.rasters;
                    rasterBins = tempRast.rasterBins;
                    % valid flag chops off from the begining and from the
                    % end numel(pspKernel)/2 bins (padLen)
                    trialSdfs = convn(rasters',pspKernel,'valid')';
                    sdfTime = rasterBins(padLen+1:end-padLen);
                    rasters = rasters(:,padLen+1:size(rasters,2)-padLen);
                    sdfMean = mean(trialSdfs,1);
                    sdfStd = std(trialSdfs,1);
                    sdfSem = sdfStd./sqrt(size(trialSdfs,1));
                    nSpikes = sum(rasters(:));
                    
                    % Gather variables for output
                    % prefix the alignedName for pivot
                    prefix = [alignedName '_'];                    
                    opts.([prefix 'alignedEvent']){roNum} = alignedEvent;
                    opts.([prefix 'alignedTimeWin']){roNum} = single(alignedTimeWin);
                    opts.([prefix 'alignTime']){roNum} = single(alignTime);
                    opts.([prefix 'firstSortByName']){roNum} = firstSortByName;
                    opts.([prefix 'firstSortByTime']){roNum} = single(firstSortByTime);
                    % may be useful to filter post-hoc
                    opts.([prefix 'nTrials'])(roNum) = numel(selTrials);
                    % may be useful to filter post-hoc
                    opts.([prefix 'nSpikes'])(roNum) = nSpikes;
                    opts.([prefix 'timeMs']){roNum} = single(sdfTime);
                    opts.([prefix 'rasters']){roNum} = rasters;
                    opts.([prefix 'sdfTsMeanStdSem']){roNum} = single([sdfTime(:) sdfMean(:) sdfStd(:) sdfSem(:)]);
                    
                end % for each alignEvent
            catch mE
                getReport(mE)
                continue
            end % try for each condition
            opts.pairedSefUnitNums{roNum}=pairedSefUnitNums;
            opts.pairedFefUnitNums{roNum}=pairedFefUnitNums;
            opts.pairedScUnitNums{roNum}=pairedScUnitNums;            
        end % for each condition
        
    end % for each unit
    fprintf(repmat('\b',1,size(s,2)));
    fprintf(' %d units done %3.2f sec\n',nUnits,toc);
end % for each area
%  optsUnitInfo table based on opts

% save spkCorrSdf.mat
fprintf('Saving output to %s...',spkCorrSdfsFile);
spkCorrSdfs = struct();
spkCorrSdfs.filterCriteria = filterCriteria;
spkCorrSdfs.unitInfo = unitInfoAll(ismember(unitInfoAll.unitNum,unique(opts.unitNum,'stable')),:);
spkCorrSdfs.SEF = opts(strcmp(opts.area,'SEF'),:);
spkCorrSdfs.FEF = opts(strcmp(opts.area,'FEF'),:);
spkCorrSdfs.SC = opts(strcmp(opts.area,'SC'),:);
save(spkCorrSdfsFile,'-v7.3','-struct','spkCorrSdfs');
fprintf(repmat('\b',1,3));
fprintf('\nDone processing in %3.2f sec\n\n',toc(stTime));

end

%%

function [pspKernel] = getPspKernelForSdf()
    % initialize the excitatory post-synaptic potential
    % Source: From Thomas Reppert: compute_spk_density_fx.m
    tau_d = 20; tau_g = 1;
    epsp = @(x) exp(-x/tau_d) .* (1 - exp(-x/tau_g));
    epsp_conv = epsp(transpose(linspace(0,199,200)));
    epsp_conv = epsp_conv * 1000/sum(epsp_conv);
    pspKernel = epsp_conv;
end


function [selTrials] = getSelectedTrials(condition,trialTypes,trialsToRemove,nTrialsThreshold)
    selTrials = trialTypes.(condition){:};
    if isempty(selTrials)
        return;
    end
    %% Mutually Exclusive trials for Choice/Timing Errors
    % For Error *Choice/*Timing ensure mutually exclusive
    otherCondition = [];
    if contains(condition,'ChoiceError')
        otherCondition = regexprep(codition,'ChoiceError','TimingError');
    elseif contains(condition, 'TimingError')
        otherCondition = regexprep(codition,'TimingError','ChoiceError');
    end
    if ~isempty(otherCondition)
        selTrials(trialTypes.(otherCondition){:}) = 0;
    end
    % First check for no selected trials: Now check for nTrials
    if isempty(selTrials) || numel(selTrials) < nTrialsThreshold
        selTrials = [];
        return;
    end
    % check if trials need to be dropped due to poor_isolation..
    if ~isempty(trialsToRemove)
        selTrials(trialsToRemove(1):trialsToRemove(2)) = 0;
    end
    % Second check for no selected trials: nTrials
    selTrials = find(selTrials);
    if numel(selTrials) < nTrialsThreshold
        selTrials = [];
    end

end

function [unitTbl] = getUniqueUnitTbl(filteredUnitsTbl)

    tempTbl = table();
    tempTbl.unitNum = [filteredUnitsTbl.X_unitNum;filteredUnitsTbl.Y_unitNum];
    tempTbl.area = [filteredUnitsTbl.X_area;filteredUnitsTbl.Y_area];
    tempTbl.sess = [filteredUnitsTbl.X_sess;filteredUnitsTbl.X_sess];
    tempTbl.rho = [filteredUnitsTbl.rho;filteredUnitsTbl.rho];
    tempTbl.rhoUnsigned = [filteredUnitsTbl.rhoUnsigned;filteredUnitsTbl.rhoUnsigned];
    tempTbl = sortrows(unique(tempTbl,'stable'),'area');
    
    uniqUnits = unique(tempTbl.unitNum, 'stable');
    
    unitTbl = table();
    for ii = 1:numel(uniqUnits)
        unitNum = uniqUnits(ii);
        idx = find(tempTbl.unitNum==unitNum);
        unitRsc = tempTbl.rho(idx);
        unitRscUnsigned = tempTbl.rhoUnsigned(idx);
        temp = table();
        
        temp.unitNum = unitNum;
        temp.sess = tempTbl.sess(idx(1));
        temp.area = tempTbl.area(idx(1));
        temp.rho = {unitRsc};
        temp.rhoUnsigned = {unitRscUnsigned};
        temp.isPlus = sum(unitRsc>0)>0;
        temp.isMinus = sum(unitRsc<0)>0;
        temp.isBoth = temp.isPlus & temp.isMinus;
        unitTbl = [unitTbl;temp]; %#ok<*AGROW>
    end
    unitTbl = sortrows(unitTbl,{'area','isBoth','isMinus','isPlus','unitNum'});
end
