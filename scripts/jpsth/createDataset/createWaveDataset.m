
wavOutputFile = 'dataProcessed/dataset/waves_SAT_DaEu.mat';
jpsthPairs = load('dataProcessed/dataset/JPSTH_PAIRS_CellInfoDB.mat');
jpsthPairs = jpsthPairs.JpsthPairCellInfoDB;

% session mat file to plexon file
sessionFiles = load('dataProcessed/dataset/SessionFiles_SAT.mat');
sessionFiles = sessionFiles.sessionFiles;
%v Need to get cell / unit nums and units
daEuPairs = innerjoin(jpsthPairs,sessionFiles,'LeftKeys',{'X_sess','X_sessNum'},...
    'RightKeys',{'session','sessNum'},...
    'LeftVariables',{'X_unitNum','X_unit','Y_unitNum','Y_unit','matDatafile'},...
    'RightVariables',{'session','sessNum','plexonFile'}...
    );

xUnits = daEuPairs(:,{'session','sessNum','X_unitNum','X_unit','matDatafile','plexonFile'});
xUnits.Properties.VariableNames = strrep(xUnits.Properties.VariableNames,'X_','');
yUnits = daEuPairs(:,{'session','sessNum','Y_unitNum','Y_unit','matDatafile','plexonFile'});
yUnits.Properties.VariableNames = strrep(yUnits.Properties.VariableNames,'Y_','');

allUnits = [xUnits;yUnits];

uniqUnits = unique(allUnits,'rows');
uniqUnits.matDatafile = regexprep(uniqUnits.matDatafile,'.*/data/','data/');
uniqUnits.unit = regexprep(uniqUnits.unit,'^(\d)([a-z])$','0$1$2');

unitsBySession = cellfun(@(x) uniqUnits(strcmp([uniqUnits.session],x),:),unique(uniqUnits.session),'UniformOutput',false);
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
save(wavOutputFile,'-mat','comment');

outChan = table();
for ii = 4:4 % 1:numel(uniqSessions)    
    sessionUnits = unitsBySession{ii};
    fprintf('Doing session %s\n',sessionUnits.session{1});
    %% Already translated data
    matFile = sessionUnits.matDatafile{1};
    matData = load(matFile,'-regexp','DSP*|TrialStart_');
    %% Plexon data and other vars needed to cut /index waveforms to trials 
    % Most recoded from RH: RH_Github_Mat_Code/Mat_Code/macTranslate
    plxFile = sessionUnits.plexonFile{1};
    plxData = readPLXFileC(plxFile,'all');
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
    
    trlEvts = arrayfun(@(x1,x2) plxEvtVal(plxEvtTs>x1 & plxEvtTs<=x2),plxTrlStartTs(1:end-1),plxTrlStartTs(2:end),'UniformOutput',false);
    trlEvts{end+1} = plxEvtVal(find(plxEvtTs>plxTrlStartTs(end)));
    trlTask = cellfun(@(x) x(find(x==taskCode,1)+1),trlEvts);
    trlHoldTime = cellfun(@(x) x(find(x==holdTimeCode,1)+1),trlEvts);
    trlTargColor = cellfun(@(x) x(find(x==targColorCode,1)+1),trlEvts);
    % from osx_breakFiles.m
    % DET_trials = find(Target_(:,9) == 1 & Target_(:,10) == 0);
    % MG_trials = find((Target_(:,9) == 1 | Target_(:,9) == 2) & Target_(:,10) > 0);
    % SEARCH_trials = find(Target_(:,9) == 2 & Target_(:,10) == 0 & Target_(:,3) < 5); %last condition ensures we were not accidentally running MG w/ a hold time of 0.  < 5 because a few days had target colors of 2 instead of 1

    detectTrls = find(trlTask == 1 & trlHoldTime == 0);
    mgTrls = find(trlTask==1 | trlTask==2 & trlHoldTime > 0);
    searchTrls = find(trlTask == 2 & trlHoldTime ==0 & trlTargColor < 5);
    
    %%
    chanNos = cellfun(@(x) str2double(x), regexp(sessionUnits.unit,'\d+','match'));    
    chanLetters = cellfun(@(x) x, regexp(sessionUnits.unit,'[a-z]$','match'));
    units = sessionUnits.unit; % 09a,10b etc
    unitNums = sessionUnits.unitNum;
    tic
    outUnits = table();
    currSess = sessionUnits.session{1};
    currSessNum = sessionUnits.sessNum(1);
    for jj = 1:numel(units)
        unitNum = unitNums(jj);
        unit = units{jj};
        chanNo = chanNos(jj);
        chanLetter = chanLetters{jj};
        outWavName = ['Unit_' num2str(unitNum,'%d')];
        outUnits.matFile{jj} = matFile;
        outUnits.plxFile{jj} = plxFile;        
        outUnits.session{jj} = currSess;
        outUnits.sessNum{jj} = currSessNum;
        outUnits.unitNum(jj) = unitNum;
        outUnits.unit{jj} = unit;
        outUnits.wavUnitName{jj} = ['WAV' num2str(chanNo,'%02d'), chanLetter];
        matUnitName = ['DSP' num2str(chanNo,'%02d'), chanLetter];
        outUnits.matUnitName{jj} = matUnitName;
        % channel data
        chData = plxData.SpikeChannels([plxData.SpikeChannels.Channel]==chanNo);
        unitsRow = strfind(letters,chanLetter)-1;
        % includes Search, Detection, Memory guided etc...
        wavTs = double(chData.Timestamps(chData.Units==unitsRow))./(plxData.ADFrequency/1000);
        wavTsIdxByTrls = arrayfun(@(x) find(wavTs>x & wavTs<=x+fixTrlLen),actualTrlStart,'UniformOutput',false);
        %wavTsIdxByTrls = arrayfun(@(x) find(wavTs>x,1):find(wavTs>x,1)+fixTrlLen,actualTrlStart,'UniformOutput',false);
        wavDat = double(chData.Waves(:,chData.Units==unitsRow))';
        wavDataByTrls = cellfun(@(x) wavDat(x,:) ,wavTsIdxByTrls,'UniformOutput',false);
        wavTsByTrls = cellfun(@(x) wavTs(x,:) ,wavTsIdxByTrls,'UniformOutput',false);
        % parse wave timestamps and waves into detection, Mg, search
        outUnits.wavDet{jj} = [];
        outUnits.wavMg{jj} = [];
        outUnits.wavSearch{jj} = [];        
        outUnits.wavDetTs{jj} = [];
        outUnits.wavMgTs{jj} = [];
        outUnits.wavSearchTs{jj} = [];
        outUnits.wavDetMean{jj} = [];
        outUnits.wavMgMean{jj} = [];
        outUnits.wavSearchMean{jj} = [];
        outUnits.wavDetStd{jj} = [];
        outUnits.wavMgStd{jj} = [];
        outUnits.wavSearchStd{jj} = [];
        
        % what is the cout of spikes by trial for the mat file and the wav file?
        outUnits.matSpkCountSearch{jj} = [];
        outUnits.wavTsCountSearch{jj} = [];
        
        if numel(detectTrls)>0
            outUnits.wavDet{jj} = arrayfun(@(x) wavDataByTrls(x),detectTrls);
            outUnits.wavDetTs{jj} = arrayfun(@(x) wavTsByTrls(x),detectTrls);
            outUnits.wavDetMean{jj} = mean(cell2mat(outUnits.wavDet{jj}));
            outUnits.wavDetStd{jj} = std(cell2mat(outUnits.wavDet{jj}));
        end
        if numel(mgTrls)>0
            outUnits.wavMg{jj} = arrayfun(@(x) wavDataByTrls(x),mgTrls);
            outUnits.wavMgTs{jj} = arrayfun(@(x) wavTsByTrls(x),mgTrls);
            outUnits.wavMgMean{jj} = mean(cell2mat(outUnits.wavMg{jj}));
            outUnits.wavMgStd{jj} = std(cell2mat(outUnits.wavMg{jj}));
        end       
        if numel(searchTrls)>0
            outUnits.wavSearch{jj} = arrayfun(@(x) wavDataByTrls(x),searchTrls);
            outUnits.wavSearchTs{jj} = arrayfun(@(x) wavTsByTrls(x),searchTrls);
            outUnits.wavSearchMean{jj} = mean(cell2mat(outUnits.wavSearch{jj}));
            outUnits.wavSearchStd{jj} = std(cell2mat(outUnits.wavSearch{jj}));
            matSpks = matData.(matUnitName);
            outUnits.matSpkCountSearch{jj} = arrayfun(@(x) sum(matSpks(x,:)>0),[1:size(matSpks,1)]');
            outUnits.wavTsCountSearch{jj} = cellfun(@(x) size(x,1),outUnits.wavSearchTs{jj});
        end
        % save every outUnit as a separate variable
        saveWavDataForUnit(wavOutputFile,['Unit_' num2str(unitNum,'%03d')],outUnits(jj,:))
    end
    toc
    %outChan = [outChan;outUnits];
    
end

function saveWavDataForUnit(oFn,unitName,unitData)
   temp.lastUnit = unitName;
   temp.(unitName)=unitData;
   save(oFn,'-append','-struct','temp');

end



