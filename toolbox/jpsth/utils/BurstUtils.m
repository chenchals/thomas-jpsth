classdef BurstUtils
    %BURSTUTILS Non-stateful burst utilities
    %   All methods are static. Computation state/results are not stored.
    %   All method calls will have to use classname prefix
    %
    
    methods (Static, Access=public)
        function outputArg = psbh(bobT, eobT, timeWin)
            %PSBH Peri Stimulus Burst Histogram
            %   Compute PSBH, PSBH-SD, PSBH-VAR for given bob, eob times.
            %   bobT : Cell array of doubles (begining of Burst Time).  The
            %      cell array must be {nTrials x 1}. Each element in the cell
            %      array is a row vector (nTimeStamps x 1) of timestamps
            %   eobT : Cell array of doubles (begining of Burst Time).  The
            %      cell array must be {nTrials x 1}. Each element in the cell
            %      array is a row vector (nTimeStamps x 1) of timestamps
            %   timeWin: [minTime maxTime] for PSBH
            %   fx:  Built-in function-handle to be used for PSBH. Valid
            %        args are: @nanmean, @nanstd, @nanvar
            
            % for each trial convert bobT_i to eobT_i
            
            %             outputArg.burstRasters = cellfun(@(x,y) burstTimes2Raster(x, y, timeWin),...
            %                 bobT, eobT, 'UniformOutput', false);
            
            outputArg.burstRasters = arrayfun(@(i) burstTimes2Raster(bobT{i}, eobT{i}, timeWin),...
                (1:size(bobT,1))', 'UniformOutput', false);
            
            % binWidth = 1 ms
            fx = @nanmean;
            outputArg.psbh = fx(cell2mat(outputArg.burstRasters));
            outputArg.rasterBins = min(timeWin):max(timeWin);
        end
        
        
        function outputArg = detectBursts(cellSpikeTimes,timeWin, varargin)
            %DETECTBURSTS
            if isempty(timeWin) % use minmax of spike times for every trial
                timeWins = cellfun(@(x) [min(x) max(x)],cellSpikeTimes,'UniformOutput',false);
                % if there are any empty timeWins, ie the trial had NaN or
                % zero spikes then, set those timeWins to [0 0]
                for ii=1:size(timeWins,1)
                    if isempty(timeWins{ii})
                        timeWins{ii,:}=[0 0];
                    end
                end
                timeWins = cell2mat(timeWins);
            elseif numel(timeWin) == 2
                timeWins = repmat(timeWin(:)',size(cellSpikeTimes,1),1);
            elseif size(cellSpikeTimes,1) == size(timeWin,1)
                timeWins = timeWin;
            elseif size(cellSpikeTimes,1) ~= size(timeWin,1)
                error(['The variable timeWin must be (a) empty or ',...
                    '(b) a vector of 2 values startT and stopT or ',...
                    '(c) a nTrials by 2 matrix']);
            end
            
            if isnumeric(timeWins)
                timeWins = mat2cell(timeWins,ones(size(timeWins,1),1),2);
            end
            
            outputArg = cellfun(@(x,y) ...
                poissBurst(x, y(1), y(2),varargin{:}),...
                cellSpikeTimes,timeWins,'UniformOutput',false);
            %fx_t = @(fn) cellfun(@(x) x.(fn),allBursts,'UniformOutput',false);
        end
        
        function outputArg = loadBursts(cellNos,fileList,blacklistedUIDs)
            outputArg = struct();
            for i = 1:numel(cellNos)
                cellNo = cellNos(i);
                cellUID =  num2str(cellNo,'UID_%04d');
                if ~sum(contains(blacklistedUIDs,cellUID))
                    fileIndex = find(~cellfun(@isempty,regexp(fileList,cellUID,'match')));
                    burstFile = fileList{fileIndex}; %#ok<FNDSB>
                    outputArg{i,1} = load(burstFile);
                end
            end
        end
                
        
        function outputArg = alignForTrials(burstTimes, varargin)
            % ALIGNFORTRIALS
            %      burstTimes: must be a cell array, if not is converted to
            %                 cell array
            
            defaultArgs = {'alignTimes', [], @isnumeric,...
                'trials', [] @isnumeric};
            % This is more complicated... what if burst spans the begining
            % or end of timeWin?
            %'timeWin', [], @(x) isnumeric(x) && numel(x)==2 && diff(x)>0,...
            
            argParser = BurstUtils.createArgParser(defaultArgs);
            if ~isempty(varargin)
                argParser.parse(varargin{:});
            end
            args = argParser.Results;
            if isnumeric(burstTimes)
                burstTimes = num2cell(burstTimes);
            end
            
            % do alignTimes first
            if numel(args.alignTimes)==1
                outputArg = cellfun(@(b) b-args.alignTimes, burstTimes,'UniformOutput', false);
            elseif numel(args.alignTimes) == size(burstTimes,1)
                alignTimes = arrayfun(@(x) {x}, args.alignTimes);
                outputArg = cellfun(@(b,t) b-t, burstTimes, alignTimes, 'UniformOutput', false);
            else
                error('Number of times in alignTimes must be 1 or equal to no of trials in burstTimes');
            end
            % if all elements of an array is  NaN replace with empty
            outputArg(cellfun(@(x) sum(isnan(x))==numel(x),outputArg))={[]};
            % Select trials, intentionally selecting AFTER alignment
            if numel(args.trials) > 0
                outputArg = outputArg(args.trials);
            end
        end
        
        % in Progress....
        function outputArg = selectBurstsInTimeWin(bobT, eobT, timeWin)
            % SELECTBURSTSINTIMEWIN Select burts in the given time window
            % for testing:
            % bt=-100:10:100;bobT=num2cell(bt(1:2:end-1));eobT=num2cell(bt(2:2:end));
            
            assert(isnumeric(timeWin) && numel(timeWin)==2 && diff(timeWin) > 0,...
                'timeWin must be a 2-element vector of doubles with diff > 0');
            
            bobInd = cellfun(@(x) find(x>=timeWin(1) & x<=timeWin(2)), bobT, 'UniformOutput', false);
            eobInd = cellfun(@(x) find(x>=timeWin(1) & x<=timeWin(2)), eobT, 'UniformOutput', false);
            % Time window includes (bobT, eobT) for all bursts
            commonInds = cellfun(@(b,e) intersect(b,e), bobInd, eobInd, 'UniformOutput', false);
            % TimWin has eobT, without bobT, first bobT is before TimeWin
            nanBobInds = find(cell2mat(cellfun(@(b,e) numel(setdiff(e,b)), bobInd, eobInd, 'UniformOutput', false)));
            % TimeWin has bobT, without eobT last eobT is after TimeWin
            nanEobInds = find(cell2mat(cellfun(@(b,e) numel(setdiff(b,e)), bobInd, eobInd, 'UniformOutput', false)));
            % get all common
            bobs = cellfun(@(b,i) b(i), bobT, bobInd,'UniformOutput',false);
            eobs = cellfun(@(e,i) e(i), eobT, eobInd,'UniformOutput',false);
            % pre-append Nan for BOB
            for j = 1:numel(nanBobInds)
                bobs{nanBobInds(j)} = [NaN bobs{nanBobInds(j)}];
            end
            % post-append Nan for EOB
            for j = 1:numel(nanEobInds)
                eobs{nanEobInds(j)} = [eobs{nanEobInds(j)} NaN];
            end
            
            outputArg.bobs = bobs;
            outputArg.eobs = eobs;
            outputArg.bothBobAndEobInds = commonInds;
            outputArg.nanBobInds = nanBobInds;
            outputArg.nanEobInds = nanEobInds;
            
        end
        
        function outputArg = convert2logical(bobT, eobT, timeWin)
            if iscell(timeWin) && size(timeWin,1) == size(bobT,1)
                outputArg = cellfun(@(b,e,t) BurstUtils.burst2logical(b,e,t), ...
                    bobT, eobT, timeWin, 'UniformOutput', false);
            else
                outputArg = cellfun(@(b,e) BurstUtils.burst2logical(b,e,timeWin), ...
                    bobT, eobT, 'UniformOutput', false);
            end
        end
                
        function saveOutput(oFile, resCellArray,varargin)
            % SAVEOUTPUT Saves output of burst analysis for a single unit
            % varargin:
            %     name-value pairs of args that are saved to file
            %var_fx = @(fn,y) cell2mat(cellfun(@(x) x.(fn),y,'UniformOutput',false))';
            var_fx = @(fn,y) cellfun(@(x) x.(fn), y,'UniformOutput',false);
            nTrials = size(resCellArray,1);
            analysisDate = datetime;
            save(oFile, 'analysisDate');
            if ~isempty(varargin)
                for i = 1:2:length(varargin)
                    o.(varargin{i}) = varargin{i+1};
                end
                save(oFile, '-append', '-struct', 'o');
            end
            
            fnames = fieldnames(resCellArray{1});
            for f = 1:numel(fnames)
                fn = fnames{f};
                switch fn
                    case {'fieldDefinitions', 'opts'}
                        t.(fn) = resCellArray{1}.(fn);
                    otherwise
                        temp = var_fx(fn,resCellArray);
                        if size(temp,2) > 1
                            % should not be the case
                            tempvals = [temp{:}];
                            t.(fn) = tempvals(:);% make a col. vector
                        else
                            t.(fn) = temp;
                        end
                end
                save(oFile,'-append','-struct', 't');
                clearvars t temp tempvals
            end
        end
        
    end
    
    methods (Static, Access=private)
        
        function argParser = createArgParser(varargin)
            argParser = inputParser();
            args = varargin{1};
            if numel(args) == 0 || mod(numel(args),3)==1
                error(['When creating argParser the no. of argumets must be greater than zero '...
                    'and must be EVEN corresponding to key-value pairs, where value is the default value']);
            end
            for  i = 1:3:numel(args)
                argParser.addParameter(args{i},args{i+1},args{i+2});
            end
            argParser.parse();
        end
        
        function [ outputArg ] = burst2logical( bobt, eobt, timeWin )
            %BURSTTIMES2RASTER Make a 0s and 1s vector for each pair of burst times
            %   Creates a zeros vector the  length of range of timeWin bins
            %   For every time-bin, add 1 if the time-bin is in burst duration  
            
            dTimeBins = min(timeWin):max(timeWin);
            outputArg = false(1,numel(dTimeBins));
            % If there are no bursts in a trial
            if isempty(bobt) || (numel(bobt)==1 && isnan(bobt))
                return
            end
            % BOB is before the timeWin
            bobt(isnan(bobt)) = timeWin(1);% replace nan with mn timeWin
            % EOB is after the timeWin
            eobt(isnan(eobt)) = timeWin(2);% replace nan with max timeWin

            burstTimes = arrayfun(@(b,e) find(dTimeBins>=b & dTimeBins<=e),bobt,eobt,'UniformOutput',false);

            burstTimes = [burstTimes{:}];
            outputArg(burstTimes) = true;
            
            
        end
        
    end
    
end

