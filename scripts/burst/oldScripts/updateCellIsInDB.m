dataDir ='/Volumes/schalllab/Users/Amir/Analysis/Mat_DataFiles';
cellInfoDbFile = 'data/Analysis/burstDB/cellInfoDB.mat';
analysisDir = 'data/Analysis/burstDB';
% Load cell inforamation database table
temp = load(cellInfoDbFile);
cellInfoTable = temp.cellInfoDB;

% update cell info database
sessions = unique(cellInfoTable.dataFile,'stable');
% get names of DSP* for each session to add to the table
temp = cellfun(@(s) who('-file',fullfile(dataDir,s),'-regexp','DSP\d+[a-h]'),...
          sessions,'UniformOutput',false);
cellIdInFile = {};
for i = 1:numel(temp)
    cellIdInFile = [cellIdInFile; temp{i}];
end
cellInfoTable.cellIdInfile=cellIdInFile;
cellInfoTable.Properties.VariableDescriptions{end} = ...
    'cellIdInfile:Ids of single inits in data file. Usually DSPXXa,..., does not include DSPxxi';
