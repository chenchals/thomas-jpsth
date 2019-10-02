% Needs ninfo_nstats_SAT.mat created by Thomas
%       and all recoded xlsx files from Rich's SAT summary excel files 
%       [monk]_SAT_colorRecode.xlsx
% Criteria for Pairs:
%      1. Only use Units the satisfy: abs(unit.visGrade) > 1 | abs(unit.moveGrade) > 1
%      2  There has to atleast 2 units in the session
% The data for all units and their criteria is in file ninfo_nstats_SAT.mat
% 
% Run createTrialTypesEventTimesDB after running this script.
% see also: CREATETRIALTYPESEVENTTIMESDB, PARSESATEXCEL

ninfo_nstat_file = 'dataProcessed/dataset/ninfo_nstats_SAT.mat';
xlsxDir = 'dataProcessed/excel';
matFilebase = '/Volumes/schalllab/data';
inRootAnalysisDir = fullfile('dataProcessed/dataset');
% if ~exist(inRootAnalysisDir,'dir')
%     mkdir(inRootAnalysisDir);
% end

monkNameMap = containers.Map({'D','E','Q','S'},{'Darwin','Euler','Quincy','Seymour'});
%% Add matlab datafile to the info table
matFileDirStruct = cellfun(@(x) dir(fullfile(matFilebase,x,'SAT/MATLAB', [x(1) '*RH_SEARCH.mat'])),...
                    monkNameMap.values,'UniformOutput',false);
matFiles = arrayfun(@(x) strcat({matFileDirStruct{x}.folder}', filesep, {matFileDirStruct{x}.name}'),...
                    1:4,'UniformOutput',false);
matFiles = vertcat(matFiles{:});

%% Add information of grid location and hemifield from recoded excel files
excelInfos = cellfun(@(x) parseSatExcel(fullfile(xlsxDir,[x '_SAT_colorRecode.xlsx'])), monkNameMap.values, 'UniformOutput', false);
excelInfos = vertcat(excelInfos{:});

%% Process info for all monks 
db = load(ninfo_nstat_file); % contains ninfo and nstats structs
cInfo = struct2table(db.ninfo);
% inner join with information from excel recoded
cInfo = innerjoin(cInfo,excelInfos,'LeftKeys',{'sess','unit'},...
                  'RightKeys',{'MatSessionName','Unit'},...
                  'RightVariables',{'MatSessionName','Hemi','Grid','SessionNotes'});

%% JPSTH summary --> how many cell pairs for each session?
JpsthPairSummary=table();
JpsthPairSummary.sess = unique(cInfo.sess);
JpsthPairSummary.nUnits = cell2mat(cellfun(@(x) sum(contains(cInfo.sess,x)),JpsthPairSummary.sess,'UniformOutput',false));
JpsthPairSummary.matFile = matFiles(contains(matFiles,JpsthPairSummary.sess));

goodVMIdx = abs(cInfo.visGrade) > 1 | abs(cInfo.moveGrade) > 1;
cellsGoodVM = cInfo(goodVMIdx,:);
cellsBySession = arrayfun(@(x) find(contains(cellsGoodVM.sess,x)), JpsthPairSummary.sess, 'UniformOutput',false);
varsForPairs = cInfo.Properties.VariableNames;
nextPairId = 0;
JpsthPairCellInfoDB = table();

for s=1:numel(cellsBySession)
    res = cellsGoodVM(cellsBySession{s},:);
    session = res.sess{1};    
    tIdx = contains(JpsthPairSummary.sess,session);
    if size(res,1) <= 1
        JpsthPairSummary.nCellsForJpsth(tIdx) = 0;
        JpsthPairSummary.nPairsJpsth(tIdx) = 0;
        continue;
    elseif size(res,1) > 1 % we have more than 1 unit
        result.CellInfoTable = cellsGoodVM(cellsBySession{s},:);
        sessName = JpsthPairSummary.sess{tIdx};
        matDatafile = JpsthPairSummary.matFile{tIdx};
        monkName = monkNameMap(char(unique(result.CellInfoTable.monkey)));
        pairRowIds = sortrows(combnk(1: size(result.CellInfoTable,1), 2),[1 2]);
        nPairs = size(pairRowIds,1);
        pairs = table();
        pairs.Pair_UID = cellstr(num2str(((1:nPairs)+ nextPairId)','PAIR_%04d'));
        
        JpsthPairSummary.nCellsForJpsth(tIdx) = size(result.CellInfoTable,1);
        JpsthPairSummary.nPairsJpsth(tIdx) = nchoosek(size(result.CellInfoTable,1),2);
        JpsthPairSummary.firstPairUID(tIdx) = pairs.Pair_UID(1);
        JpsthPairSummary.lastPairUID(tIdx) = pairs.Pair_UID(end);
        
        for v = 1:numel(varsForPairs)
            cName = varsForPairs{v};
            pairs.(['X_' cName]) = result.CellInfoTable.(cName)(pairRowIds(:,1));
            pairs.(['Y_' cName]) = result.CellInfoTable.(cName)(pairRowIds(:,2));
        end
        pairs.matDatafile = repmat({matDatafile},nPairs,1);          
        nextPairId = nextPairId + nPairs;
        JpsthPairCellInfoDB = [JpsthPairCellInfoDB;pairs]; %#ok<AGROW>
        tempSumm =table();
        tempSumm.sessionName = sessName;
    end
    result.PairInfoTable = pairs;
end
save(fullfile(inRootAnalysisDir,'JPSTH_PAIRS_CellInfoDB.mat'),'JpsthPairCellInfoDB');
save(fullfile(inRootAnalysisDir,'JPSTH-PAIR-Summary.mat'), 'JpsthPairSummary');
writetable(JpsthPairSummary,fullfile(inRootAnalysisDir,'JPSTH-PAIR-Summary.csv'))


