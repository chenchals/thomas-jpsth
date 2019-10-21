function [] = updateWaveDataset(waveformDir)
%UPDATEWAVEDATASET Summary of this function goes here

wavFiles = dir(fullfile(waveformDir,'Unit_*'));
wavFiles = strcat({wavFiles.folder}','/',{wavFiles.name}');
for ii = 1:numel(wavFiles)
    wavFile = wavFiles{ii};
    wavTbl = load(wavFile);
    fns = fieldnames(wavTbl);
    unitName = char(fns(contains(fns,'Unit_')));
    wavNames = {'wavDet', 'wavSearch', 'wavMg'};
    wavWidthNames = {'wavWidthDet', 'wavWidthSearch', 'wavWidthMg'};
    
    for jj = 1:numel(wavNames)
        wavTbl.(unitName).(wavWidthNames{jj}) = getWaveformWidths(wavTbl.(unitName).(wavNames{jj}));
    end
    updateWavDataForUnit(wavFile,wavTbl);
    
end

end

function updateWavDataForUnit(oFn,wavTbl)
   temp.lastUpdated=datestr(datetime());
   fns = fieldnames(wavTbl);
   for ii = 1:numel(fns)
    temp.(fns{ii}) = wavTbl.(fns{ii});
   end
   fprintf('Updating file %s\n',oFn);
   save(oFn,'-struct','temp');
   %save(oFn,'-append','-struct','temp');
end