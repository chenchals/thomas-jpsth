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

% Modifications:
% 10/03/2019 : Removed the criteria for filtering for only vis/mov units

useVizMovUnits = false;
ninfo_nstat_file = 'dataProcessed/dataset/ninfo_nstats_SAT.mat';
xlsxDir = 'dataProcessed/excel';
matAndPlxFiles = 'dataProcessed/dataset/SessionFiles_SAT.mat';
inRootAnalysisDir = fullfile('dataProcessed/dataset');
% if ~exist(inRootAnalysisDir,'dir')
%     mkdir(inRootAnalysisDir);
% end

monkNameMap = containers.Map({'D','E','Q','S'},{'Darwin','Euler','Quincy','Seymour'});
%% Add matlab datafile and plx datafile  to the info table
sessMatPlxFiles = load(matAndPlxFiles);
sessMatPlxFiles = sessMatPlxFiles.sessionFiles;
%% Add information of grid location and hemifield from recoded excel files
excelInfos = cellfun(@(x) parseSatExcel(fullfile(xlsxDir,[x '_SAT_colorRecode.xlsx'])), monkNameMap.values, 'UniformOutput', false);
excelInfos = vertcat(excelInfos{:});

%% Process info for all monks 
db = load(ninfo_nstat_file); % contains ninfo and nstats structs
cInfo = struct2table(db.ninfo);

% convert grid location from ap-ml to 
% +(plus) for anterior -(minus) for posterior
% +(plus) for medial -(minus) for lateral
fx_GridConvert = @(gloc) fliplr(eval([ '[' regexprep(fliplr(gloc),{'a','p','m','l'},{'+','-','+','-'}) ']']));
gridLoc = regexprep(lower(excelInfos.Grid),{'-','nan',' '},{''});
excelInfos.GridAP_ML = cellfun(@(x) fx_GridConvert(x),gridLoc,'UniformOutput',false);

newDepth = str2double(excelInfos.Depth);
d0Idx = isnan(str2double(excelInfos.Depth)) & ~isnan(str2double(excelInfos.Depth0));
newDepth(d0Idx) = str2double(excelInfos.Depth0(d0Idx));
excelInfos.newDepth = newDepth;

% inner join with information from excel recoded
cInfo = innerjoin(cInfo,excelInfos,'LeftKeys',{'sess','unit'},...
                  'RightKeys',{'MatSessionName','Unit'},...
                  'RightVariables',{'MatSessionName','Hemi','Grid','GridAP_ML',...
                  'Depth','Depth0','newDepth','SessionNotes'});

%% JPSTH summary --> how many cell pairs for each session?
JpsthPairSummary=table();
JpsthPairSummary.sess = unique(cInfo.sess);
JpsthPairSummary.nUnits = cell2mat(cellfun(@(x) sum(contains(cInfo.sess,x)),JpsthPairSummary.sess,'UniformOutput',false));
temp = innerjoin(sessMatPlxFiles,JpsthPairSummary);
JpsthPairSummary.matDatafile = temp.matDatafile;
JpsthPairSummary.plxDatafile = temp.plxDatafile;

if useVizMovUnits
    goodVMIdx = abs(cInfo.visGrade) > 1 | abs(cInfo.moveGrade) > 1; %#ok<UNRCH>
    goodUnits = cInfo(goodVMIdx,:);
else
    goodUnits = cInfo;
end

cellsBySession = arrayfun(@(x) find(contains(goodUnits.sess,x)), JpsthPairSummary.sess, 'UniformOutput',false);
varsForPairs = cInfo.Properties.VariableNames;
nextPairId = 0;
JpsthPairCellInfoDB = table();

%% Ensure the following is true for area1 vs area2
% for cross area pairs: always have first SEF on x-axis, then FEF on
% x-axis, then SC on x-axis; NSEFN is always on y-axis
pairXYarea = {
          {'SEF' 'SEF'}
          {'SEF' 'FEF'}
          {'SEF' 'SC'}
          {'SEF' 'NSEFN'}
          {'FEF' 'FEF'}
          {'FEF' 'SC'}
          {'FEF' 'NSEFN'}
          {'SC' 'SC'}
          {'SC' 'NSEFN'}         
          {'NSEFN' 'NSEFN'}
          };
% concated strings
pairXYarea = cellfun(@(x) [x{:}],pairXYarea,'UniformOutput',false);
          
for s=1:numel(cellsBySession)
    res = goodUnits(cellsBySession{s},:);
    session = res.sess{1};    
    tIdx = contains(JpsthPairSummary.sess,session);
    if size(res,1) <= 1
        JpsthPairSummary.nCellsForJpsth(tIdx) = 0;
        JpsthPairSummary.nPairsJpsth(tIdx) = 0;
        continue;
    elseif size(res,1) > 1 % we have more than 1 unit
        result.CellInfoTable = goodUnits(cellsBySession{s},:);
        sessName = JpsthPairSummary.sess{tIdx};
        matDatafile = JpsthPairSummary.matDatafile{tIdx};
        plxDatafile = JpsthPairSummary.plxDatafile{tIdx};        
        
        monkName = monkNameMap(char(unique(result.CellInfoTable.monkey)));
        pairRowIds = sortrows(combnk(1: size(result.CellInfoTable,1), 2),[1 2]);
        nPairs = size(pairRowIds,1);
        pairs = table();
        pairs.Pair_UID = cellstr(num2str(((1:nPairs)+ nextPairId)','PAIR_%04d'));
        
        JpsthPairSummary.nCellsForJpsth(tIdx) = size(result.CellInfoTable,1);
        JpsthPairSummary.nPairsJpsth(tIdx) = nchoosek(size(result.CellInfoTable,1),2);
        JpsthPairSummary.firstPairUID(tIdx) = pairs.Pair_UID(1);
        JpsthPairSummary.lastPairUID(tIdx) = pairs.Pair_UID(end);
        
        % add flag to check if we need to swap rowIds so as to get X-Unit
        % and Y-Uint too be congruent with pairXYAreas defined above
        % the pairRowIds(:,3) is the swar flag, swap col1 and 2 for that
        % row
        pairRowIds(:,3) = arrayfun(@(x) ...
            sum(strcmp(pairXYarea,[result.CellInfoTable.area{pairRowIds(x,:)}]))== 0,...
            1:nPairs)';
        pairRowIds(:,4:5) = pairRowIds(:,1:2);
        swapCols = find(pairRowIds(:,3));
        for zz = 1:numel(swapCols)
            pairRowIds(swapCols(zz),1:2) = fliplr(pairRowIds(swapCols(zz),1:2));
        end
        for v = 1:numel(varsForPairs)
            cName = varsForPairs{v};
            pairs.(['X_' cName]) = result.CellInfoTable.(cName)(pairRowIds(:,1));
            pairs.(['Y_' cName]) = result.CellInfoTable.(cName)(pairRowIds(:,2));
        end
        pairs.XY_Dist = arrayfun(@(x)...
            getElectrodeDistance(pairs.X_GridAP_ML{x},pairs.Y_GridAP_ML{x},...
                                pairs.X_newDepth(x),pairs.Y_newDepth(x)),...
                                (1:size(pairs,1))','UniformOutput',false);
        pairs.isOnSameChannel = arrayfun(@(x) ...
            isequal(regexp(pairs.X_unit{x},'(\d{1,2})','tokens'),...
            regexp(pairs.Y_unit{x},'(\d{1,2})','tokens')), (1:size(pairs))');
        
        pairs.matDatafile = repmat({matDatafile},nPairs,1);          
        pairs.plxDatafile = repmat({plxDatafile},nPairs,1);          
        JpsthPairCellInfoDB = [JpsthPairCellInfoDB;pairs]; %#ok<AGROW>
        tempSumm =table();
        tempSumm.sessionName = sessName;
        nextPairId = nextPairId + nPairs;
    end
    result.PairInfoTable = pairs;
end
save(fullfile(inRootAnalysisDir,'JPSTH_PAIRS_CellInfoDB.mat'),'JpsthPairCellInfoDB');
save(fullfile(inRootAnalysisDir,'JPSTH-PAIR-Summary.mat'), 'JpsthPairSummary');
writetable(JpsthPairSummary,fullfile(inRootAnalysisDir,'JPSTH-PAIR-Summary.csv'))


function [eDist] = getElectrodeDistance(xGridLoc,yGridLoc,xDepth,yDepth)
% grid locs must be numeric [AP, ML]
% [+1 = 1a, -1 = 1p] [+1 = 1m, -1 = 1L]
% depth in microns --> convert to mm

    if isempty(xGridLoc) || isempty(yGridLoc) || sum(isnan([xDepth yDepth]))
        eDist = NaN;
    elseif isequal(xGridLoc,yGridLoc)
        eDist = 0;
    else
        pointA = [xGridLoc(:); xDepth/1000];
        pointB = [yGridLoc(:); yDepth/1000];
        eDist =  sqrt(sum((pointA-pointB).^2));
    end

end
