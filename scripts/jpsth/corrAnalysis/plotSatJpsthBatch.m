
jpsthDataBaseDir = 'dataProcessed/analysis/JPSTH';

pairDirs = dir(fullfile(jpsthDataBaseDir,'jpsth_*'));
pairDirs = strcat({pairDirs.folder}',filesep,{pairDirs.name}');
tic
saveFigFlag = 1;
for pd = 1:numel(pairDirs)
    pairDir = pairDirs{pd};
    pdfOutputDir = fullfile(pairDir,'pdf');
    if ~exist(pdfOutputDir,'dir')
      mkdir(pdfOutputDir);
    end
    
    jpsthPairFiles = dir(fullfile(pairDir,'JPSTH-PAIR_*.mat'));
    jpsthPairFiles = strcat({jpsthPairFiles.folder}',filesep,{jpsthPairFiles.name}');
    % cant do parfor as plot gets clipped
    for pf = 1:numel(jpsthPairFiles)
        jpsthPairFile = jpsthPairFiles{pf};
        [~,temp] = fileparts(jpsthPairFile);
        fprintf('Plotting JPSTH data for file : %s ...',temp);
        plotSatJpsth(jpsthPairFile,pdfOutputDir,saveFigFlag);
        fprintf('done!\n');
    end
end
toc
