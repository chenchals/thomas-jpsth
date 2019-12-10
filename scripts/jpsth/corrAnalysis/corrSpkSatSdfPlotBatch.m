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

availEpochToUse = {'Baseline','PostSaccade'};
availOutcomes = {'Correct','ErrorChoice','ErrorTiming'};
areas = {'SEF','FEF','SC'};
baseSatSdfFile = 'dataProcessed/analysis/11-18-2019/spkCorr/summary/spkCorrSdfs';
baseSatSdfPdfDir = 'dataProcessed/analysis/11-18-2019/spkCorr/summary/sdfPlots';

for ae = 1:numel(availEpochToUse)
    epochToUse = availEpochToUse{ae};
    
    for ao = 1:numel(availOutcomes)
        outcome = availOutcomes{ao};
         
        basePdfDir = fullfile(baseSatSdfPdfDir,sprintf('using_%s_%s',outcome,epochToUse));
        if ~exist(basePdfDir,'dir')
            mkdir(basePdfDir);
        end
        
        satSDfFiles = dir(sprintf('%s_%s_%s.mat',baseSatSdfFile,outcome,epochToUse));
        satSdfFiles = strcat({satSDfFiles.folder},filesep,{satSDfFiles.name})';
        sdfsTbl = table();
        unitInfosTbl = table();
        
        currOutcome = outcome;
        currFile = satSdfFiles{contains(satSdfFiles,currOutcome)};
        temp = load(currFile);
        tempUnitInfo = temp.unitInfo;
        tempUnitInfo.Properties.RowNames = {};
        unitInfosTbl = [unitInfosTbl;tempUnitInfo]; %#ok<*AGROW>
        for a = 1:numel(areas)
            currArea = areas{a};
            aT = temp.(currArea);
            aT.filterCriteria = repmat({temp.filterCriteria},size(aT,1),1);
            sdfsTbl = [sdfsTbl;aT];
        end % for each area
        
        uniqUnits = unique(sdfsTbl.unitNum);
        tic
        for u = 5:numel(uniqUnits)
            tic
            % call function: corrSpkSatSdfPlot(unitSdfsTbl,unitInfoTbl,pdfFilename)
            currUnitNum = uniqUnits(u);
            fprintf('Doing unit num %03d...\n', currUnitNum);
            currUnitSdfsTbl = sdfsTbl(sdfsTbl.unitNum==currUnitNum,:);
            currUnitInfoTbl = unitInfosTbl(unitInfosTbl.unitNum==currUnitNum,:);
            currUnitInfoTbl = currUnitInfoTbl(1,:);
            currPdfFilename = sprintf('%sUnit_%03d_%s.pdf',[basePdfDir filesep],currUnitNum,currUnitInfoTbl.area{1});
            corrSpkSatSdfPlot(currUnitSdfsTbl,currUnitInfoTbl,currPdfFilename);
            fprintf('Unit num %03d [%d of %d] units done %3.2f sec\n',currUnitNum, u, numel(uniqUnits),toc);
            
        end
    end
end


