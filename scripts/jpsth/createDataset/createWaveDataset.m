
wavOutputDir = 'dataProcessed/dataset/waves2';
if ~exist(wavOutputDir,'dir')
    mkdir(wavOutputDir);
end
% load ninfo
ninfo = load('dataProcessed/dataset/ninfo_nstats_SAT.mat');
ninfo = struct2table(ninfo.ninfo);

% session mat file to plexon file
sessionFiles = load('dataProcessed/dataset/SessionFiles_SAT.mat');
sessionFiles = sessionFiles.sessionFiles;

ninfo = innerjoin(ninfo,sessionFiles,'LeftKeys',{'sess','sessNum'},...
    'RightKeys',{'sess','sessNum'});
ninfo.unit = cellfun(@(x) regexprep(x,'^(\d)([a-z])$','0$1$2'),ninfo.unit,'UniformOutput',false);

uniqSessions = unique(ninfo.sess);
unitsBySession = cellfun(@(x) ninfo(strcmp([ninfo.sess],x),:),uniqSessions,'UniformOutput',false);
% map cluster to row no. in plexon file
% row 0='i',row 1='a'...
letters = 'iabcdefghjklmnopqrstuvwxyz';

%from RH: RH_Github_Mat_Code/Mat_Code/macTranslate
%The Task code follows the strobe "3003" and will be:
%0: fixation    1: detection    2: search/MG
taskCode = 3003;
taskCodedDef = {{0 'fixation'};{1 'detection'};{2 'search'}};
% Hold time: to check of MG or Search: value follw stribe "3021"
holdTimeCode = 3021;
trlStartCode = 1666;
targColorCode = 3008;
% EventChannel. name = Strobed
%strobeIdx = 17; % channel 257
pdIdx = 2; % channel 2 photodiode
baselineTime = 3500;%mseconds
fixTrlLen = 6000; % mqseconds

%% Create an empty output file to addend wave data for each unit
comment = {'Waveform data for Darwin and Euler SAT data, (not all units)'... 
          'See ninfo_nstats_SAT.mat (ninfo) for unit details'};
%save(wavOutputFile,'-mat','comment');
warning('off')
tic
for ii = 1:numel(uniqSessions)
    sessionUnits = unitsBySession{ii};
    fprintf('Doing session [%d of %d] %s\n',ii,numel(uniqSessions),sessionUnits.sess{1});
    %% Already translated data
    matFile = sessionUnits.matDatafile{1};
    fprintf('Loading matFile ')
    matData = load(matFile,'-regexp','DSP*|TrialStart_');
    %% Plexon data and other vars needed to cut /index waveforms to trials 
    % Most recoded from RH: RH_Github_Mat_Code/Mat_Code/macTranslate
    plxFile = sessionUnits.plxDatafile{1};
    fprintf('plexon file...')
    plxData = readPLXFileC(plxFile,'all');
    fprintf('Done!\n');
    strobeIdx = strcmp({plxData.EventChannels.Name},'Strobed');
    plxEvtTs = double(plxData.EventChannels(strobeIdx).Timestamps)./plxData.ADFrequency;
    plxEvtVal = plxData.EventChannels(strobeIdx).Values;
    plxPdTs = double(plxData.EventChannels(pdIdx).Timestamps)./plxData.ADFrequency;
    isPd = numel(plxPdTs)>0;
    plxTrlStartTs = plxEvtTs(plxEvtVal==trlStartCode);
    targOnPd = plxTrlStartTs;
    if isPd
        targOnPd = arrayfun(@(x) plxPdTs(find(plxPdTs>x,1)), plxTrlStartTs);
    end
    actualTrlStart = floor(targOnPd*1000) - baselineTime;    
    trlEvts = arrayfun(@(x1,x2) plxEvtVal(plxEvtTs>x1 & plxEvtTs<x2-1),plxTrlStartTs(1:end-1),plxTrlStartTs(2:end),'UniformOutput',false);
    %trlEvts{end+1} = plxEvtVal(find(plxEvtTs>plxTrlStartTs(end)));
    trlEvts{end+1} = plxEvtVal(plxEvtTs>plxTrlStartTs(end)); %#ok<SAGROW>
    % for those trialStarts, where the task has not started.. usually the
    % last trial, when experiment stops at the momemt TrialStart_ code is
    % sent....Then we will not find index for task code "3003"
    trlTaskIdx = cellfun(@(x) find(x==taskCode,1)+1,trlEvts,'UniformOutput',false);
    trlHoldTimeIdx = cellfun(@(x) find(x==holdTimeCode,1)+1,trlEvts,'UniformOutput',false);
    % find valid trials
    validTrlIdx = find(~cellfun(@(x,y) isempty(x)|isempty(y),trlTaskIdx,trlHoldTimeIdx));
    % prune timestamps and events
    actualTrlStart = actualTrlStart(validTrlIdx);
    trlEvts = trlEvts(validTrlIdx);
    trlTaskIdx = trlTaskIdx(validTrlIdx);
    trlTask = cellfun(@(x,y) x(y),trlEvts,trlTaskIdx);
    trlHoldTime = cellfun(@(x) x(find(x==holdTimeCode,1)+1),trlEvts);
    trlTargColor = cellfun(@(x) x(find(x==targColorCode,1)+1),trlEvts);
    % from osx_breakFiles.m
    % DET_trials = find(Target_(:,9) == 1 & Target_(:,10) == 0);
    % MG_trials = find((Target_(:,9) == 1 | Target_(:,9) == 2) & Target_(:,10) > 0);
    % SEARCH_trials = find(Target_(:,9) == 2 & Target_(:,10) == 0 & Target_(:,3) < 5); 
    %      last condition ensures we were not accidentally running MG w/ a hold time of 0. 
    %       < 5 because a few days had target colors of 2 instead of 1
    detectTrls = find(trlTask == 1 & trlHoldTime == 0);
    mgTrls = find(trlTask==1 | trlTask==2 & trlHoldTime > 0);
    searchTrls = find(trlTask == 2 & trlHoldTime ==0 & trlTargColor < 5);
    
    %%
    chanNos = cellfun(@(x) str2double(x), regexp(sessionUnits.unit,'\d+','match'));    
    chanLetters = cellfun(@(x) x, regexp(sessionUnits.unit,'[a-z]$','match'));
    units = sessionUnits.unit; % 09a,10b etc
    unitNums = sessionUnits.unitNum;
    tic
    currSess = sessionUnits.sess{1};
    currSessNum = sessionUnits.sessNum(1);
    parfor jj = 1:numel(units)
        outUnits = table();
        unitNum = unitNums(jj);
        unit = units{jj};
        %fprintf('......Unit %s [%d of %d]\n',unit,jj,numel(units)); %#ok<PFBNS>
        chanNo = chanNos(jj);
        chanLetter = chanLetters{jj};
        outWavName = ['Unit_' num2str(unitNum,'%d')];
        outUnits.matFile{1} = matFile;
        outUnits.plxFile{1} = plxFile;        
        outUnits.session{1} = currSess;
        outUnits.sessNum{1} = currSessNum;
        outUnits.unitNum{1} = unitNum;
        outUnits.unit{1} = unit;
        outUnits.wavUnitName{1} = ['WAV' num2str(chanNo,'%02d'), chanLetter];
        matUnitName = ['DSP' num2str(chanNo,'%02d'), chanLetter];
        outUnits.matUnitName{1} = matUnitName;
        % channel data
        chData = plxData.SpikeChannels([plxData.SpikeChannels.Channel]==chanNo);
        unitsRow = strfind(letters,chanLetter)-1;
        % includes Search, Detection, Memory guided etc...
        wavTs = double(chData.Timestamps(chData.Units==unitsRow))./(plxData.ADFrequency/1000);
        wavTsIdxByTrls = arrayfun(@(x) find(wavTs>x & wavTs<=x+fixTrlLen),actualTrlStart,'UniformOutput',false);
        %wavTsIdxByTrls = arrayfun(@(x) find(wavTs>x,1):find(wavTs>x,1)+fixTrlLen,actualTrlStart,'UniformOutput',false);
        wavDat = single(chData.Waves(:,chData.Units==unitsRow))';
        wavDataByTrls = cellfun(@(x) wavDat(x,:) ,wavTsIdxByTrls,'UniformOutput',false);
        wavTsByTrls = cellfun(@(x) wavTs(x,:) ,wavTsIdxByTrls,'UniformOutput',false);
        % parse wave timestamps and waves into detection, Mg, search
        outUnits.wavDet{1} = [];
        outUnits.wavMg{1} = [];
        outUnits.wavSearch{1} = [];        
        outUnits.wavDetTs{1} = [];
        outUnits.wavMgTs{1} = [];
        outUnits.wavSearchTs{1} = [];
        outUnits.wavDetMean{1} = [];
        outUnits.wavMgMean{1} = [];
        outUnits.wavSearchMean{1} = [];
        outUnits.wavDetStd{1} = [];
        outUnits.wavMgStd{1} = [];
        outUnits.wavSearchStd{1} = [];
        
        % what is the cout of spikes by trial for the mat file and the wav file?
        outUnits.matSpkCountSearch{1} = [];
        outUnits.wavTsCountSearch{1} = [];
        outUnits.isEqualSpkCountsSearchMatWav{1} = [];
        outUnits.diffSpkCountsSearchMinMax{1} = [];
        
        if numel(detectTrls)>0
            outUnits.wavDet{1} = arrayfun(@(x) wavDataByTrls(x),detectTrls);
            outUnits.wavDetTs{1} = arrayfun(@(x) wavTsByTrls(x),detectTrls);
            outUnits.wavDetMean{1} = mean(cell2mat(outUnits.wavDet{1}));
            outUnits.wavDetStd{1} = std(cell2mat(outUnits.wavDet{1}));
        end
        if numel(mgTrls)>0
            outUnits.wavMg{1} = arrayfun(@(x) wavDataByTrls(x),mgTrls);
            outUnits.wavMgTs{1} = arrayfun(@(x) wavTsByTrls(x),mgTrls);
            outUnits.wavMgMean{1} = mean(cell2mat(outUnits.wavMg{1}));
            outUnits.wavMgStd{1} = std(cell2mat(outUnits.wavMg{1}));
        end       
        if numel(searchTrls)>0
            outUnits.wavSearch{1} = arrayfun(@(x) wavDataByTrls(x),searchTrls);
            % these are trial timestamps from trial start
            outUnits.wavSearchTs{1} = arrayfun(@(x) round(wavTsByTrls{x}-actualTrlStart(x)),searchTrls,'UniformOutput',false);
            outUnits.wavSearchMean{1} = mean(cell2mat(outUnits.wavSearch{1}));
            outUnits.wavSearchStd{1} = std(cell2mat(outUnits.wavSearch{1}));
            matSpks = matData.(matUnitName);
            outUnits.matSpkCountSearch{1} = arrayfun(@(x) sum(matSpks(x,:)>0),(1:size(matSpks,1))');
            outUnits.wavTsCountSearch{1} = cellfun(@(x) size(x,1),outUnits.wavSearchTs{1});
            outUnits.isEqualSpkCountsSearchMatWav{1} = isequal(outUnits.matSpkCountSearch{1},outUnits.wavTsCountSearch{1});
            d = outUnits.matSpkCountSearch{1} - outUnits.wavTsCountSearch{1};
            outUnits.diffSpkCountsSearchMinMax{1} = minmax((outUnits.matSpkCountSearch{1} - outUnits.wavTsCountSearch{1})');
        end
        % save every outUnit as a separate variable
        oUnitName = ['Unit_' num2str(unitNum,'%03d')];
        saveWavDataForUnit(fullfile(wavOutputDir,oUnitName),oUnitName,outUnits)
    end
    toc  
end
toc


function saveWavDataForUnit(oFn,unitName,unitData)
   temp.lastSaved=datestr(datetime());
   temp.(unitName)=unitData;
   save(oFn,'-struct','temp');
   %save(oFn,'-append','-struct','temp');
end



