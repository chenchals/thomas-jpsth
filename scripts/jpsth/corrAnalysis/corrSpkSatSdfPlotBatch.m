% Create a unique unit table for units chosen for given epoch that show
% significant correlations and plot SDFs for each unique unit. 
% Dataset(s) already created by running corrSpkSatSdf.m, that creates the
% following files for different epochs:
% Baseline epoch: 
%   summary/spkCorrSdfs_Correct_Baseline.mat
%   summary/spkCorrSdfs_ErrorChoice_Baseline.mat
%   summary/spkCorrSdfs_ErrorTiming_Baseline.mat
% PostSaccade epoch: 
%   summary/spkCorrSdfs_Correct_PostSaccade.mat
%   summary/spkCorrSdfs_ErrorChoice_PostSaccade.mat
%   summary/spkCorrSdfs_ErrorTiming_PostSaccade.mat
% see also CORRSPKSATSDF

epochToUse = 'Baseline';
outcomes = {'Correct','ErrorChoice','ErrorTiming'};
areas = {'SEF','FEF','SC'};
baseSatSdfFile = 'dataProcessed/analysis/11-18-2019/spkCorr/summary/spkCorrSdfs_';
baseSatSdfPdfDir = 'dataProcessed/analysis/11-18-2019/spkCorr/summary/sdfPlots';
basePdfDir = fullfile(baseSatSdfPdfDir,['using_',epochToUse]);
if ~exist(basePdfDir,'dir')
    mkdir(basePdfDir);
end

satSDfFiles = dir([baseSatSdfFile '*_' epochToUse '.mat']);
satSdfFiles = strcat({satSDfFiles.folder},filesep,{satSDfFiles.name})';
sdfsTbl = table();
unitInfosTbl = table();
for o = 1:numel(outcomes)
    currOutcome = outcomes{o};
    currFile = satSdfFiles{contains(satSdfFiles,currOutcome)};
    temp = load(currFile);
    tempUnitInfo = temp.unitInfo;
    tempUnitInfo.Properties.RowNames = {};
    unitInfosTbl = [unitInfosTbl;tempUnitInfo]; %#ok<*AGROW>
    allAreasTbl = table();    
    for a = 1:numel(areas)
        currArea = areas{a};
        aT = temp.(currArea);
        aT.filterCriteria = repmat({temp.filterCriteria},size(aT,1),1);
        sdfsTbl = [sdfsTbl;aT];         
    end % for each area
end % for each outcome

uniqUnits = unique(sdfsTbl.unitNum);

for u = 1:numel(uniqUnits)
    % call function: corrSpkSatSdfPlot(unitSdfsTbl,unitInfoTbl,pdfFilename)
    currUnitNum = uniqUnits(u);
    currUnitSdfsTbl = sdfsTbl(sdfsTbl.unitNum==currUnitNum,:);
    currUnitInfoTbl = unitInfosTbl(unitInfosTbl.unitNum==currUnitNum,:);
    currUnitInfoTbl = currUnitInfoTbl(1,:);
    pdfFilename = sprintf('%sUnit_%03d_%s.pdf',[basePdfDir filesep],currUnitNum,currUnitInfoTbl.area{1});
    
end



