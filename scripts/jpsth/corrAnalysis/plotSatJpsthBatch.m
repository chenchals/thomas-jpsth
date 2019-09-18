
jpsthDataBaseDir = 'dataProcessed/analysis/JPSTH';

pairDirs = {
    'jpsth_SEF-FEF/mat'
    'jpsth_SEF-SC/mat'
    'jpsth_FEF-SC/mat'
    'jpsth_SEF-SEF/mat'
    'jpsth_FEF-FEF/mat'
    'jpsth_SC-SC/mat'
    };

pairDirs = strcat(jpsthDataBaseDir,filesep,pairDirs);

tic
saveFigFlag = 1;
for pd = 1:numel(pairDirs)
    pairDir = pairDirs{pd};
    pdfOutputDir = fullfile(pairDir,'pdf2');
    if ~exist(pdfOutputDir,'dir')
      mkdir(pdfOutputDir);
    end
    
    jpsthPairFiles = dir(fullfile(pairDir,'JPSTH-PAIR_*.mat'));
    jpsthPairFiles = strcat({jpsthPairFiles.folder}',filesep,{jpsthPairFiles.name}');
    % cant do parfor as plot gets clipped
    for pf = 1:numel(jpsthPairFiles)
        try
        jpsthPairFile = jpsthPairFiles{pf};
        [~,temp] = fileparts(jpsthPairFile);
  
            fprintf('Plotting JPSTH data for file : %s ...',temp);
            plotSatJpsth(jpsthPairFile,pdfOutputDir,saveFigFlag);

        fprintf('done!\n');
        catch me
            fprintf('****** ERROR*****\n')
            disp(me)
            continue
        end
    end
end
toc
