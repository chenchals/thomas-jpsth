function [oStruct]=poissBurst(inputTrain, inStartT, inStopT,varargin)
% POISSBURST Modified from p_burst.m (from Hanes, Thompson, and Schall 1995). See
% http://www.psy.vanderbilt.edu/faculty/schall/pdfs/reofpres.pdf
% http://www.psy.vanderbilt.edu/faculty/schall/matlab/P_burst.zip)
%Returns Spike Index of BOB, EOB and Surprise Of Burst
%The spike train is a series of timestamps
%Function p_burst
%NOTE: Calls POISSCDF function in the STATISTICS toolbox
%INPUTS:
%		InTrain: Timestamps of the spike train.
%		inStartT: Window start time
%       inStopT: window stop time
%         If no varargs, then spike count between inStartT and inStopT is
%         used to calculate average firing rate MU for complete spike train.
%         Burst results may not be valid outside limits. 
%       varargin: name-value pairs of variable arguments
%         averageFr: Mu to be used for analysis
%         useWindowForMu: boolean, 
%                         true [default]: use spike count between inStartT and
%                                         inStopT for computing Mu 
%                         false: use all spikes in inTrain for computing Mu
%         plotBursts: boolean, 
%                         false [default]: do not plot.  
%                         true: plot bursts and raster. uses existing
%                               figure. Only meaningfule for debugging.                        
%         minSpikesInBurst: scalar 2 [default]. Min number spikes in burst  
%         jitterSpikeTimes: boolean, useful when there are simultaneous
%                                    spike times
%                           false [default] do not jitter spike times
%                           true add jitter time to spike times
%         anchorTime: 
%         maxExtraSpikes:
%         maxExtraTime:
%OUTPUTS: outStruct
%     fieldDefinitions: definitions of each field as below
%     spkT: spike times or train
%     spkTWin: spike times in time window
%     timeWin: time window used to compute mu for burst
%     bobT: begining of burst time
%     eobT: end of burst time
%     dobT: duration of burst
%     nsdb: number of spikes during burst
%     ibiT: inter burst interval time
%     nsdibi: number of spikes during inter burst interval
%     frdb: firing rate during burst
%     frdbAvg: firing rate during burst average
%     frdibi: firing rate during inter burst interval
%     frdibiAvg: firing rate during inter burst interval average
%     frdTWin: firing rate during time window
%     frSpkT: firing rate for spike train input
%PROBLEMS:
%		If you get an error: One or more output arguments
%		not assigned during call to 'distchck'
%		Change the spike Train dimension by passing
%		a transpose(Spike Train)
%CALLING:
%[oStruct]=p_burst(InTrain, StartT, StopT, varargin)
%
% First Version Jan 5, 1998: S.Chenchal Rao
% Modifications Mar 9, 2018: S.Chenchal Rao
%    use relative spike times
%    use argument parser for input control of parameters
%
%
%******************User Parameters************************
opts = parseInputArgs(varargin);

%%%%%%%%%%%%%%%% Function %%%%%%%%%%%%%%%%%%%%%%%%
tempInTrain = inputTrain;
MaxXT=opts.maxExtraTime; % 30;%Max Xtra Time
MaxXS=opts.maxExtraSpikes;% 10;%Max Xtra Spikes
MinSPInBurst=opts.minSpikesInBurst;% 2;%Minimum spkes in a Burst
Anchor=opts.anchorTime; %50;%Anchor Time
Significance=opts.Significance;
opts.UserSI=-log(Significance);
UserSI=opts.UserSI;
Tol=1e-300;
if opts.jitterSpikeTimes
    jitterForSimultaneousSpikes = 1e-20;
    [InTrain, jitterTimes] = jitterDuplicateTimestamps(inputTrain,jitterForSimultaneousSpikes);
else
    jitterTimes = 0;
    InTrain = inputTrain;
end
%****************Spike Train Properties*****************
if(size(InTrain,1))>1
    InTrain=InTrain';
end
%%%%%%%%%%%%%%% Modifications - 2018 %%%%%%%%%%%%%%%%%%%
if(~inStopT),inStopT=max(InTrain);end
if(~inStartT),inStartT=min(InTrain);end
% Ensure start and stop times are not negative
% Offset window time if negative
if inStartT > inStopT
    error('Start time MUST be lower than stop time');
end

offsetSpkTime = min(min(InTrain),min(inStartT));

StartT = inStartT - offsetSpkTime;
StopT = inStopT - offsetSpkTime;
% Add offestSpkTime to InTrain, such that 0 corresponds to StartT
SPT = InTrain - offsetSpkTime;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use spike times >
useSpkIndexes = SPT >0 & SPT < Inf;
SPT = SPT(useSpkIndexes)';
if numel(jitterTimes) > numel(SPT)
  jitterTimes = jitterTimes(useSpkIndexes);
end

if(~isempty(SPT))
    StopT=max(0,StopT);%traps NaN
    StartT=max(0,StartT);%traps NaN
    % if(~StopT),StopT=max(SPT);end
    % if(~StartT),StartT=min(SPT);end
    minT=min(StartT,StopT);
    maxT=max(StartT,StopT);
    StartT=minT;StopT=maxT;
    Duration=StopT-StartT;
    %Average Spike Rate MU
    if opts.averageFr > 0
        % since mu given in Hz, divide by 1000
        MU = opts.averageFr/1000;
    elseif opts.useWindowForMu
        MU = (length(find(SPT>= StartT & SPT <= StopT)))/(Duration);
    elseif ~opts.useWindowForMu
        MU = (length(SPT)-1)/(max(SPT)-min(SPT));
    end
    %if there are no spikes in the spike train or
    %n spikes in train is less than 4
    if(MU==0 || length(SPT)<=4)
        BOB=[];
        EOB=[];
        SOB=[];
        oStruct = cleanOutput(false);
        return
    end
else
    BOB=[];
    EOB=[];
    SOB=[];
    oStruct = cleanOutput(false);
    return
end
MaxSpikes=length(SPT);
ISI=diff(SPT);%Inter Spike Intervals
ISI = ISI + Tol;
%################Parameter Initializations##############
%******Flags, Counters, Indices and Rel. OPs************
ISBURST=0;
BNo=1;CurrBNo=0;%Burst No for indexing
MinBOB=5000;%Dummy
MaxEOB=0;%Dummy
OldBOB=0;OldEOB=0;
%PrevBOB=0;PrevEOB=0;
CurrBOB=0;CurrEOB=0;
%CurrEOB_XS=0;CurrXS=0;CurrXT=0;%Vars for Extra Time
%MaxSI=0;
%*******************Iteration vars**********************
%Iterate=1;
IC=1;%Iteration Counter
%FromI=0;ToI=0; %Current portion of Spike Train
FspAB=1;%First spike After Burst
%Temp=0;
Done=0;
%******************Output Arguments*********************
%BOB=single([]);EOB=single([]);SOB=single([]);%MUST be set to EMPTY

BOB = zeros(1,5000,'single');
EOB = zeros(1,5000,'single');
SOB = zeros(1,5000,'single');
if opts.useGpu
 gpuDevice(); % reset gpu device
 fprintf('Using GPU\n');
 SPT = gpuArray(SPT);
 ISI = gpuArray(ISI);
 BOB = gpuArray(zeros(1,5000,'single'));
 EOB = gpuArray(zeros(1,5000,'single'));
 SOB = gpuArray(zeros(1,5000,'single'));
end
%########################################################
while(FspAB <= MaxSpikes-1 || ~Done)
    Iterate=1;
    while(Iterate)
        %***************FIND EOB****************************
        if (IC==1)
            FromI=FspAB;
        else
            FromI=CurrBOB;
        end
        cISI=cumsum(ISI(FromI:length(ISI)))+Anchor;
        Prob=(poisscdf((1:length(cISI))', cISI.*MU))+Tol;
        Prob=(1-Prob);SI=-log(Prob);
        %find index that maximizes SI
        Temp=(find(diff(SI)<0));
        if (length(Temp)>1)
            %         fprintf('Temp: \n')
            %         disp(Temp)
            CurrEOB=Temp(min(find(Temp >1)));
            if(CurrEOB==1 && IC==1)
                FspAB=Temp(min(find(diff(Temp) >MinSPInBurst)))+1;
                if (isempty(FspAB) || FspAB>=MaxSpikes-MinSPInBurst)
                    Done=1;
                    ISBURST=0;
                    Iterate=0;
                end
                break %Break out of Iteration
            end
            CurrEOB_XS=0;Curr_XT=0;
            %check for extra spikes and extra time
            %Find index to next MaxSI
            CurrEOB_XS=min(find(SI > SI(CurrEOB)));
            if(~isempty(CurrEOB_XS))
                CurrXS=CurrEOB_XS-CurrEOB-1;
                CurrXT=SPT(CurrEOB_XS)-SPT(CurrEOB);
                if(CurrXS <= MaxXS && CurrXT <= MaxXT)
                    CurrEOB=CurrEOB_XS;
                    CurrEOB_XS=0;
                end
            end
        elseif (length(Temp)==1 && ~isempty(Temp))
            CurrEOB=Temp;
        elseif(isempty(Temp))
            %Index to the max
            CurrEOB=length(SI);
        else
        end
        CurrEOB=CurrEOB+FromI;
        
        %***************FIND BOB****************************
        ToI=FspAB;        
        BSPT=SPT(CurrEOB:-1:ToI);
        %         fprintf('BSPT: \n')
        %         disp(BSPT)
        if numel(BSPT==1)
            BSPT = [BSPT; BSPT+Tol];
        end
        
        cISI=cumsum(abs(diff(BSPT)))+Anchor;
        Prob=poisscdf([1:length(cISI)]',cISI.*MU)+Tol;
        Prob=(1-Prob);
        SI=-log(Prob);
        Temp=max(find(SI==max(SI)));
        %What if Temp = [];??
        if (isempty(Temp))
            disp('ERROR finding BOB ==> Temp is EMPTY')
            [Temp]
        end
        CurrBOB=find(SPT==BSPT(Temp+1));
        % Modify for simultaniety
        CurrBOB = min(CurrBOB);
        
        %if diff between CurrEOB and recently found BOB is < criteria
        if (CurrEOB-CurrBOB) < MinSPInBurst
            ISBURST=0;
            Iterate=0;
            FspAB=CurrEOB+1;
            FromI=FspAB;
            break;%Break out of Iteration loop
        end
        %Check for convergence among 3 values for each
        if(OldBOB==CurrBOB && OldEOB==CurrEOB )
            Iterate=0;
            ISBURST=1;
            break
        else
            IC=IC+1;
            PrevBOB=OldBOB;
            PrevEOB=OldEOB;
            OldBOB=CurrBOB;
            OldEOB=CurrEOB;
        end
        %find max and min
        CurrBNo=max(CurrBNo,BNo);
        if(CurrBNo)
            MaxEOB=max(CurrEOB,MaxEOB);
            MinBOB=min(CurrBOB,MinBOB);
        end
        if(IC==10)
            MaxEOB=CurrEOB;MinBOB=CurrBOB;
        elseif(IC>20)
            CurrEOB=MaxEOB;CurrBOB=MinBOB;
            Iterate=0;
            ISBURST=1;
        else
        end
        %disp(['End of Iteration cycle : ',num2str(IC-1)])
        
    end%while(Iterate)
    if(ISBURST)
        %CHECK CRITERIA
        MaxSI=-log(1-(poisscdf(CurrEOB-CurrBOB,(SPT(CurrEOB)-SPT(CurrBOB))*MU)));
        if(CurrEOB-CurrBOB >= MinSPInBurst && MaxSI>UserSI)
            BOB(BNo)=CurrBOB;
            EOB(BNo)=CurrEOB;
            SOB(BNo)=MaxSI;
            BNo=BNo+1;
        end
        FspAB=CurrEOB+1;
        CurrBOB=FspAB;
        IC=1;
    else
        IC=1;
    end
    %Check limits
    if(FspAB >= MaxSpikes-2)
        Done=1;
        Iterate=0;
        break
    end
end%while(FspAB <= MaxSpikes-1)

% Subtract offsetSpkTime to get back original times for SPT
if opts.useGpu
  SPT = gather(SPT);
  BOB = gather(BOB);
  BOB = gather(BOB);
  SOB = gather(SOB);
  gpuDevice(); % clear gpuDevice
end
SPT = SPT + offsetSpkTime - jitterTimes;
BOBT = SPT(BOB(BOB>0))';
EOBT = SPT(EOB(EOB>0))';
SOB = SOB(SOB>0);

% only spkTimes in window
SPTWin = SPT(SPT>=inStartT & SPT<=inStopT);
if opts.plotBursts
    plotIt;
end
% Duration of burst
DOBT = EOBT - BOBT;
% Inter burst interval Time
IBIT = abs(EOBT(1:end-1)-BOBT(2:end));
% Firing rate during burst (FRDB)
% Below not work if SPT are sub millisecond, so
% FRDB = 1000*((EOB-BOB)+1)./(EOBT-BOBT) ;
% Count *all* spikes during burst * 1000
FRDB = 1000*(arrayfun(@(b,e) numel(SPT(SPT>=b & SPT<=e)),BOBT,EOBT))./(EOBT-BOBT);
% number of spikes during burst
NSDB = arrayfun(@(b,e) numel(SPT(SPT>=b & SPT<=e)),BOBT,EOBT);
% Firing rate Inter Burst
% Count *all* spikes during INTER burst * 1000
FRDIBI =  1000*(arrayfun(@(b,e) numel(SPT(SPT>=b & SPT<=e)),EOBT(1:end-1),BOBT(2:end)))./IBIT;
% number of spikes Inter burst
NSDIBI = arrayfun(@(b,e) numel(SPT(SPT>=b & SPT<=e)),EOBT(1:end-1),BOBT(2:end));


oStruct = cleanOutput(true);

    function oStruct = cleanOutput(analyzed)
        
        % Output as struct?
        oStruct.fieldDefinitions = {
            % 'spkT:spike times or train';
            'spkTWin:spike times in time window';
            'opts:options used for poissBurst call';
            'timeWin:time window used to compute mu for burst';
            'nb:number of bursts';
            'bobT:begining of burst time';
            'eobT:end of burst time';
            'dobT:duration of burst';
            'sob:surprise of burst';
            'nsdb:number of spikes during burst';
            'ibiT:inter burst interval time';
            'nsdibi:number of spikes during inter burst interval';
            'frdb:firing rate during burst';
            'frdbAvg:firing rate during burst average';
            'frdibi:firing rate during inter burst interval';
            'frdibiAvg:firing rate during inter burst interval average';
            'frdTWin:firing rate during time window';
            'frSpkT:firing rate for spike train input';
            };
        oStruct.opts = opts;
        oStruct.timeWin = [inStartT inStopT];
        if ~analyzed
            % oStruct.spkT = [];
            oStruct.spkTWin = [];
            oStruct.nb = 0;
            oStruct.bobT = [];
            oStruct.eobT = [];
            oStruct.dobT = [];
            oStruct.sob = [];
            oStruct.nsdb = [];
            oStruct.ibiT = [];
            oStruct.nsdibi = [];
            oStruct.frdb = [];
            oStruct.frdbAvg = [];
            oStruct.frdibi = [];
            oStruct.frdibiAvg = [];
            oStruct.frdTWin = [];
            oStruct.frSpkT = [];
            if exist('SPT','var')
                %oStruct.spkT = SPT;
                oStruct.frSpkT = NaN;
                if ~isempty(SPT)
                    oStruct.frSpkT = 1000*numel(SPT)/range(SPT);
                end
            end
            if exist('SPTWin','var')
                %oStruct.spkTWin = SPTWin;
                oStruct.frdTWin = NaN;
                if ~isempty(SPT)
                    oStruct.frdTWin = 1000*numel(SPTWin)/range([inStartT,inStopT]);
                end
            end
        else
            %oStruct.spkT = SPT;
            oStruct.spkTWin = SPTWin;
            oStruct.nb = numel(BOBT);
            oStruct.bobT = BOBT;
            oStruct.eobT = EOBT;
            oStruct.dobT = DOBT;
            oStruct.sob = SOB;
            oStruct.nsdb = NSDB;
            oStruct.ibiT = IBIT;
            oStruct.nsdibi = NSDIBI;
            oStruct.frdb = FRDB;
            oStruct.frdbAvg = nanmean(FRDB);
            oStruct.frdibi = FRDIBI;
            oStruct.frdibiAvg = nanmean(FRDIBI);
            oStruct.frdTWin = 1000*numel(SPTWin)/range([inStartT,inStopT]);
            oStruct.frSpkT = 1000*numel(SPT)/range(SPT);
        end
        % replace all []/Inf fields with NaN
        fields = fieldnames(oStruct);
        for f = 1:numel(fields)
            val = oStruct.(fields{f});
            if isnumeric(val) && (isempty(val) || sum(isinf(val)))
                oStruct.(fields{f}) = NaN;
            end
        end
    end

    function [] = plotIt()
        %****************************DISPLAY***********************************
        clf
        % do not create a new figure, else there will be 1 fig per call, it
        % will be too slow and also crash!
        %figure
        colors_='rbm';
        colors_=repmat(colors_,1,ceil(length(EOBT)/length(colors_)));
        %It takes less time to put '+' than draw line/ticks
        %h=plot(SPT,1,'k+');
        set(line([SPT SPT],[0.95 1.05]),'color',[0 0 0],'linewidth',0.1);
        hold on
        %axis([min(SPT) max(SPT) 0 2])
        axis([inStartT inStopT 0 2])
        set(findobj('type','axes'),'color',[1 1 1],'ytick',[],'box','on')
        for burstIndex=1:length(EOBT)
            plot(SPT(SPT>=BOBT(burstIndex) & SPT<=EOBT(burstIndex)),1.1,strcat(colors_(burstIndex),'+'));
        end
        drawnow
    end
end

function [ args ] = parseInputArgs(varargin)
   
   argParser = inputParser();
   argParser.addParameter('jitterSpikeTimes', false); % default
   argParser.addParameter('useWindowForMu', true) % default
   argParser.addParameter('minSpikesInBurst', 2); % default
   argParser.addParameter('averageFr', 0); % default
   argParser.addParameter('maxExtraTime', 30); % default
   argParser.addParameter('maxExtraSpikes', 10); % default
   argParser.addParameter('anchorTime', 50); % default
   argParser.addParameter('plotBursts', false); % default  
   argParser.addParameter('useGpu', false); % default  
   argParser.addParameter('Significance', 0.05); % default  
   argParser.parse();
   if ~isempty(varargin{1})
      argParser.parse(varargin{1}{:});
   end
   args = argParser.Results;
end



function [oTrain, jitterTimes] = jitterDuplicateTimestamps(inTrain, jitter)
    % add increments of jitter (1e-6) msecs to spikes that are simultaneous
    % ex. spktimes [200 200 200 201 203] ==>
    % after jitter = [200 200.000001 200.000002 201 203]

    jitterTimes = randperm(length(inTrain))'*jitter;
    oTrain = sort(inTrain + jitterTimes);
end


% Call for SpikeTimes.saccade
% allCells = arrayfun(@(c) SpikeTimes.saccade(:,c),1:29,'UniformOutput',false);
% OStruct = arrayfun(@(x) poissBurst(x{1},-200,800),c1,'UniformOutput',false);
% allCells = cellfun(@(c) c, SpikeTimes.saccade,'UniformOutput',false);

% fx_allCells = @(allC) cellfun(@(c) c, allC,'UniformOutput',false);
% OStructAllCells = arrayfun(@(x) poissBurst(x{1},-200,800),fx_allCells(SpikeTimes.saccade),'UniformOutput',false);
% allOStruct = arrayfun(@(c) arrayfun(@(x) poissBurst(x{1},-200,800),c,'UniformOutput',false),SpikeTimes.saccade,'UniformOutput',false);

% to get to allBurst fields:
% cellIndex = 5;
%allBobT = arrayfun(@(x) res(1).allBursts{x,cellIndex}.bobT, (1:size(allBursts,1))','UniformOutput',false);
