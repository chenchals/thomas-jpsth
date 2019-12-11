function [outTbl] = getUnitSatSdf(useUnit,evntTimes,useTrials,useAlignment)
%GETUNITSATSDF Summary of this function goes here
%   Detailed explanation goes here

%     useUnit.unitNum = unitNum;
%     useUnit.spkTimes = spikesSat{unitNum}';
%     evntTimes = sessionEventTimes(strcmp(sessionEventTimes.session,sess),:);
%     useTrials.trialTypes = sessionTrialTypes(strcmp(sessionTrialTypes.session,sess),:);
%     useTrials.trialsToRemove = unitInfoAll.trRemSAT{unitInfoAll.unitNum==unitNum};
%     useTrials.minTrialCount = 1;
%     useAlignment.events = {'CueOn','SaccadePrimary','RewardTime'};
%     useAlignment.timeWins = {[-600 400],[-200 600],[-100 400]};
%     useAlignment.names = {'Visual','PostSaccade','PostReward'};
%     useAlignment.sortEventNames = {'SaccadePrimary','SaccadePrimary','SaccadePrimary'};
%

unitNum = useUnit.unitNum;
spkTimes = useUnit.spkTimes; %spikesSat{unitNum}';

trialTypes = useTrials.trialTypes;%  sessionTrialTypes(strcmp(sessionTrialTypes.session,sess),:);
trialsToRemove = useTrials.trialsToRemove;% unitInfoAll.trRemSAT{unitInfoAll.unitNum==unitNum};
minTrialCount = useTrials.minTrialCount;

alignEvents = useAlignment.events; %{'CueOn','SaccadePrimary','RewardTime'};
alignTimeWins = useAlignment.timeWins; % {[-600 400],[-200 600],[-100 400]};
alignNames = useAlignment.names; %{'Visual','PostSaccade','PostReward'};
sortEventNames = useAlignment.sortEventNames; %{'SaccadePrimary','SaccadePrimary','SaccadePrimary'};


%conditions
satConditions = {
    'AccurateCorrect';'AccurateErrorChoice';'AccurateErrorTiming';
    'FastCorrect';    'FastErrorChoice';    'FastErrorTiming'
    };

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

%%
roNum = 0;
outTbl = table();
for cc = 1:numel(satConditions)
    % SDF for each condition,
    % some condition may not exist...so use try...catch
    try
        condition = satConditions{cc};
        
        % Get trials for this unit discounting the trials to remove
        selTrials = getSelectedTrials(condition,trialTypes,trialsToRemove,minTrialCount);
        if isempty(selTrials)
            continue;
        end
        % Increment row counter for output - pivot alignEvents
        % There is 1 row per condition. The column names capture
        % ALIGNED_NAME_[raters,sdf etc]
        roNum = roNum + 1;
        outTbl.unitNum(roNum) = unitNum;
        outTbl.condition{roNum} = condition;
        % pivot all aligned events. That is make separate columns
        % for each aligned event by prefixing the columnName with
        % aligned name
        for evId = 1:numel(alignEvents)
            alignedEvent = alignEvents{evId};
            alignedTimeWin = alignTimeWins{evId}; %#ok<*PFBNS>
            alignedTimeWinPad = alignedTimeWin + padTimeWin;
            alignedName = alignNames{evId};
            % sort trials
            firstSortByName = sortEventNames{evId};
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
            alignedSpkTimes = SpikeUtils.alignSpikeTimes(spkTimes(selTrials),alignTime, alignedTimeWinPad);
            tempRast = SpikeUtils.rasters(alignedSpkTimes,alignedTimeWinPad);
            rasters = tempRast.rasters;
            rasterBins = tempRast.rasterBins;
            trialSdfs = convn(rasters',pspKernel,'full')';
            trialSdfs = trialSdfs(:,padLen:range(alignedTimeWin)+padLen);
            sdfTime = rasterBins(padLen+1:end-padLen);
            rasters = rasters(:,padLen+1:size(rasters,2)-padLen);
            sdfMean = mean(trialSdfs,1);
            sdfStd = std(trialSdfs,1);
            sdfSem = sdfStd./sqrt(size(trialSdfs,1));
            nSpikes = sum(rasters(:));
            
            % Gather variables for output
            % prefix the alignedName for pivot
            prefix = [alignedName '_'];
            outTbl.([prefix 'alignedEvent']){roNum} = alignedEvent;
            outTbl.([prefix 'alignedTimeWin']){roNum} = single(alignedTimeWin);
            outTbl.([prefix 'alignTime']){roNum} = single(alignTime);
            outTbl.([prefix 'firstSortByName']){roNum} = firstSortByName;
            outTbl.([prefix 'firstSortByTime']){roNum} = single(firstSortByTime);
            % may be useful to filter post-hoc
            outTbl.([prefix 'nTrials'])(roNum) = numel(selTrials);
            % may be useful to filter post-hoc
            outTbl.([prefix 'nSpikes'])(roNum) = nSpikes;
            outTbl.([prefix 'timeMs']){roNum} = single(sdfTime);
            outTbl.([prefix 'rasters']){roNum} = rasters;
            outTbl.([prefix 'sdfTsMeanStdSem']){roNum} = single([sdfTime(:) sdfMean(:) sdfStd(:) sdfSem(:)]);
            
        end % for each alignEvent
    catch mE
        getReport(mE)
        continue
    end % try for each condition
end % for each condition




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


function [pspKernel] = getPspKernelForSdf()
    % initialize the excitatory post-synaptic potential
    % Source: From Thomas Reppert: compute_spk_density_fx.m
    tau_d = 20; tau_g = 1;
    epsp = @(x) exp(-x/tau_d) .* (1 - exp(-x/tau_g));
    epsp_conv = epsp(transpose(linspace(0,199,200)));
    epsp_conv = epsp_conv * 1000/sum(epsp_conv);
    pspKernel = epsp_conv;
end

