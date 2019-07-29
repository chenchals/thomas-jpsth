
    allDaData ='/mnt/teba/Users/Thomas/0-chenchal/Info/Darwin.mat'; % not useful
    analysisDir = '/mnt/teba/Users/Thomas/0-chenchal/BurstAnalysis/burstAlignedDBDarwin';
    monk = 'Da'; % used as suffix for the info
    useTrialTypes = {'SAT' 'MG'};
    
    %% Processing Logic
    if ~exist(analysisDir,'dir')
        mkdir(analysisDir);
    end
    
    temp = load(allDaData);
    %% Load cell information and add specific fields
    CellInfoDB = temp.(['ninfo' monk]);% DET,MG,SAT
    nCells = numel(CellInfoDB);% number of cells
    CellInfoDB = struct2table(CellInfoDB); %
    t.UID = num2str((1:nCells)','UID_%04d');
    t.datafile = CellInfoDB.sesh;
    t.monk = arrayfun(@(x) monk,(1:nCells)','UniformOutput',false);
    t.SessionNo = CellInfoDB.snum;
    t.cellIdInFile = strcat('DSP',CellInfoDB.unit);
    CellInfoDB = [struct2table(t) CellInfoDB];
    save(fullfile(analysisDir,'CellInfoDB.mat'),'CellInfoDB');
    
    
    % Load Spike Time data 
    spikeData = temp.(['spikes' monk]);%DET,MG,SAT
    % fieldnames are the different trial types DET, MG, SAT
    f = fieldnames(spikeData);
    spikeData = struct2table(spikeData);
    spikeData.Properties.VariableNames = f;
    % clear temp var
    clearvars temp t
    
    cellNos = 1:nCells;
     % for each cell -> ensure spike times is a celll array of {nTrials x 1 }
     
    parfor c = 1:nCells
         for tt = 1:numel(useTrialTypes) 
             fprintf('Detecting bursts for Cell #%d and trialType %s \n',c,useTrialTypes{tt});
             trialType = useTrialTypes{tt};
             cellInfo = CellInfoDB(c,:);
             sessionNo = cellInfo.SessionNo;
             % Save 1 file for each trial type
             fileName = char(join({cellInfo.UID,'session',num2str(sessionNo,'%03d'),'aligned', trialType},'_'));
             analysisFile = fullfile(analysisDir, [fileName '.mat']);
             % Convert to colum of trials with row
             spkTimes = spikeData.(trialType){c}';
             timeWins = cellfun(@(x) minmax(x),spkTimes,'UniformOutput',false);
             oBursts = BurstUtils.detectBursts(spkTimes,timeWins);
             BurstUtils.saveOutput(analysisFile,oBursts,'cellInfo',cellInfo);
             waitForFileOutput(analysisFile);
             fprintf('wrote file %s\n\n',analysisFile);
             % since we are using aligned timestamps need ot create a
             % logical array for isBursting
             createIsBurstingAndSave(analysisFile);
             fprintf('updated file %s with isBursting logical\n\n',analysisFile);
         end
     end
     
     
     function  createIsBurstingAndSave(analysisFile)
         oBursts = load(analysisFile);
         oBursts.isBursting = BurstUtils.convert2logical(oBursts.bobT,oBursts.eobT,[-2500 3500]);
         %overwrite
         save(analysisFile, '-struct', 'oBursts');
     end
     
     function waitForFileOutput(analysisFile)
          while(~exist(analysisFile,'file'))
              pause(0.2);
          end
     end

     
