% The code in createSatJpsthDataset.m has area matching issues for SEF
% tyhe specific lines of offending code are:
%% Filter cell pairs for the areas of interest
% jpsthCellPairs = jpsthCellPairs(...
%     ((contains(jpsthCellPairs.X_area,area1) & contains(jpsthCellPairs.Y_area,area2)) ...
%    | (contains(jpsthCellPairs.X_area,area2) & contains(jpsthCellPairs.Y_area,area1))),...
%      :);
% The function contains should not be used as for example to match:
%   SEF vs SEF --> contains will also match area NSEFN 
%                  since NSEFN contains SEF
% So anywhere we match for SEF, we will also match for NSEFN which results directory
%  SEF-SEF to also contain SEF-NSEF, as well as NEFN-NSEFN files
%  SEF-FEF to also contain NSEF-FEF
%  SEF-SC  to also contain NSEF-SC
%  SEF_NSEF to also contain NSEF-NSEF
% 
% To fix this we delete the offending files from the respective directories

areaPairs ={
    'SEF','SEF'
    'SEF','FEF'
    'SEF','SC'
    'SEF','NSEFN'
    'FEF','FEF' % should be OK
    'FEF','SC' % should be OK
    'FEF','NSEFN' % should be OK
    'SC','SC' % should be OK
    'SC','NSEFN' % should be OK
    'NSEFN','NSEFN' % should be OK
    };
% sanity check:
JpsthCellPairs = load('dataProcessed/dataset/JPSTH_PAIRS_CellInfoDB.mat');
jpsthCellPairs = JpsthCellPairs.JpsthPairCellInfoDB;

fx_incorrectPairs = @(area1,area2) ...
    ((contains(jpsthCellPairs.X_area,area1) & contains(jpsthCellPairs.Y_area,area2)) ...
    | (contains(jpsthCellPairs.X_area,area2) & contains(jpsthCellPairs.Y_area,area1)));
fx_correctPairs = @(area1,area2) ...
    ((strcmp(jpsthCellPairs.X_area,area1) & strcmp(jpsthCellPairs.Y_area,area2)) ...
    | (strcmp(jpsthCellPairs.X_area,area2) & strcmp(jpsthCellPairs.Y_area,area1)));
jpsthBaseDir = 'dataProcessed/analysis/JPSTH-10ms';
spkCorrBaseDir = 'dataProcessed/analysis/spkCorr';
incorrPairsTbl = {};
corrPairsTbl = {};
incorrP = 0;
corrP = 0;
for p = 1:size(areaPairs,1)
    area1 = areaPairs{p,1};
    area2 = areaPairs{p,2};
    areaDir = [area1 '-' area2];
    jpsthDelDir = fullfile(jpsthBaseDir, ['jpsth_' areaDir],'deleted');
    jpsthDirMat = fullfile(jpsthBaseDir, ['jpsth_' areaDir],'mat');
    jpsthDirPdf = fullfile(jpsthBaseDir, ['jpsth_' areaDir],'pdf');
    spkCorrDelDir = fullfile(spkCorrBaseDir, ['spkCorr_' areaDir],'deleted');
    spkCorrDirMat = fullfile(spkCorrBaseDir, ['spkCorr_' areaDir],'mat');
    spkCorrDirPdf = fullfile(spkCorrBaseDir, ['spkCorr_' areaDir],'pdf');
    % get pairs numbers
    incorrectPairs = find(fx_incorrectPairs(area1,area2));
    correctPairs = find(fx_correctPairs(area1,area2));
    % 
    matFiles = dir(fullfile(jpsthDirMat,'*.mat'));
    pdfFiles = dir(fullfile(jpsthDirPdf,'*.pdf'));
    if isequal(incorrectPairs,correctPairs)
        corrP = corrP + 1;
        fprintf('areaPair %s-%s is OK.  Incorrect %d, Correct %d\n',area1,area2,numel(incorrectPairs),numel(correctPairs));
         temp = sortrows(jpsthCellPairs(correctPairs,{'Pair_UID','X_area','Y_area'}),{'X_area','Y_area'});
         temp.areaDir = repmat(areaDir,size(temp,1),1);
         corrPairsTbl{corrP,1} = temp;
    else
         incorrP = incorrP + 1;
         fprintf('areaPair %s-%s is NOT OK\n',area1,area2);
         fprintf('\t\tIncorrect Pairs %d\n',numel(incorrectPairs));
         fprintf('\t\tCorrect Pairs %d\n',numel(correctPairs));
         setDiffIncorrCorr = setdiff(incorrectPairs,correctPairs);
         setDiffCorrIncorr = setdiff(incorrectPairs,correctPairs);
         fprintf('\t\tSet diff (incorrect,correct) %d\n',numel(setDiffIncorrCorr));
         fprintf('\t\tSet diff (correct,incorrect) %d\n',numel(setDiffCorrIncorr));
         temp = sortrows(jpsthCellPairs(setDiffIncorrCorr,{'Pair_UID','X_area','Y_area'}),{'X_area','Y_area'});
         temp.areaDir = repmat(areaDir,size(temp,1),1);
         incorrPairsTbl{incorrP,1} = temp;
%          
%          if ~exist(jpsthDelDir,'dir')
%              mkdir(jpsthDelDir);
%          end
%          if ~exist(spkCorrDelDir,'dir')
%              mkdir(spkCorrDelDir);
%          end
%          move files from mat and pdf dirs to deleted dir       
%          for jj = 1:size(temp,1)
%              pairId = temp.Pair_UID{jj};
%              for jpsth
%                 src = [jpsthDirMat,'/*',pairId,'*']; 
%                 dest = jpsthDelDir;
%                 try
%                 movefile(src,dest);
%                 catch me
%                     disp(me)
%                 end
%                 src = [jpsthDirPdf,'/*',pairId,'*'];
%                 try
%                 movefile(src,dest);
%                 catch me
%                     disp(me)
%                 end
%                 for spkCorr
%                 src = [spkCorrDirMat,'/*',pairId,'*']; 
%                 dest = spkCorrDelDir;
%                 try
%                 movefile(src,dest);
%                 catch me
%                     disp(me)
%                 end
%                 src = [spkCorrDirPdf,'/*',pairId,'*']; 
%                 try
%                 movefile(src,dest);
%                 catch me
%                     disp(me)
%                 end            
%          end % move end
     end
    
end


