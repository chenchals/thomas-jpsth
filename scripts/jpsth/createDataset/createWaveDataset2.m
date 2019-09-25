
% read waveform from plexon file for every spike time of a given
% session-unit 
dataDir = '/Users/subravcr/teba';
wavOutputDir = 'dataProcessedLocal/dataset/waves2';
if ~exist(wavOutputDir,'dir')
    mkdir(wavOutputDir);
end

% get neuron info
ninfo = load('dataProcessed/dataset/ninfo_nstats_SAT.mat','ninfo');
ninfo = struct2table(ninfo.ninfo);

% session files for plexon...
sessionFiles = load('dataProcessed/dataset/SessionFiles_SAT.mat');
sessionFiles = sessionFiles.sessionFiles;

unitsBySession = cellfun(@(x) ninfo(strcmp(ninfo.sess,x),:), sessionFiles.session,'UniformOutput',false);

for s = 1:numel(unitsBySession)
    units = unitsBySession{s};
    chanNos = cellfun(@(x) str2double(x), regexp(units.unit,'\d+','match'));    
    chanLetters = cellfun(@(x) x, regexp(units.unit,'[a-z]$','match'));
    sess = units.sess{1};
    sessFiles = sessionFiles(strcmp(sessionFiles.session,sess),:);
    plxFile = fullfile(dataDir,sessFiles.plexonFile{1});
    matFileBase = regexprep(plxFile,{'/Plexon' '/Sorted' '/Unsorted' '\.plx'},{'/Matlab' '' '' ''});    
    % load the plexon file...
    plx = readPLXFileC(plxFile,'all');
    fsKHz = plx.ADFrequency/1000;
    % load plexon data
    plx=struct2table(plx.SpikeChannels);
     
    % load translated data from mat files
    mat.DET = [];
    if exist([matFileBase '_DET.mat'],'file')
        mat.DET = load([matFileBase '_DET.mat'],'-regexp','^DSP*|TrialStart_');
    end
    
    mat.MG = [];
    if exist([matFileBase '_MG.mat'],'file')
        mat.MG = load([matFileBase '_MG.mat'],'-regexp','^DSP*|TrialStart_');
    end
    
    mat.SEARCH = [];
    if exist([matFileBase '_SEARCH.mat'],'file')
        mat.SEARCH = load([matFileBase '_SEARCH.mat'],'-regexp','^DSP*|TrialStart_');
    end
    
    for uu = 1:size(units,1)
        unit = units(uu,:);
        chanNo = chanNos(uu);
        chanLet = chanLetters{uu};
        unitName = ['DSP' num2str(chanNo,'%02d') chanLet];
        u = table();
        u.session{1} = sess;
        u.sessNum(1) = unit.sessNum;       
        u.unitNum(1) = unit.unitNum;
        u.unit{1} = unit.unit{1};
        % load a unit and TrialStart_ for a session
        unitSpks = mat.SEARCH.(unitName);
        trlStrt = mat.SEARCH.TrialStart_(:,1);
        unitSpks(unitSpks==0) = NaN;
        unitSpks = unitSpks + trlStrt;
        unitSpks = arrayfun(@(x) unitSpks(x,~isnan(unitSpks(x,:))),(1:size(unitSpks,1)),'UniformOutput',false);
        
        wavDat = single(plx.Waves{chanNo})';
        wavTs = single(plx.Timestamps{chanNo})./fsKHz;
        wavSpks = arrayfun(@(x) wavTs(wavTs>x & wavTs<=x+6000)',trlStrt,'UniformOutput',false);       
        wavDat = arrayfun(@(x) wavDat(wavTs>x & wavTs<=x+6000,:),trlStrt,'UniformOutput',false);
        
        wIdx = cell(numel(wavSpks),1);
        waves = cell(numel(wavSpks),1);
        for t = 1: numel(wavSpks)
            wv = wavDat{t};
            w = wavSpks{t};
            u = unitSpks{t};
            wIdx{t,1} = cell2mat(arrayfun(@(x) find(abs(round(w-x))<=1,1),u,'UniformOutput',false));
            %wv = w(wIdx{t})';
            waves{t} = wv(wIdx{t},:);
        end
        
    end
end
