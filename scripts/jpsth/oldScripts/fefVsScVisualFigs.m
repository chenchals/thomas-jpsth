
% Do for all fefvsScVisual filesfiles...
jpDir = '/Volumes/schalllab/Users/Chenchal/JPSTH/FEF_SC_Visual_1ms';
d=dir(fullfile(jpDir,'PAIR*'));
dataFiles = strcat({d.folder}',filesep,{d.name}');

for jj=1:numel(dataFiles)
    datafile=dataFiles{jj};
     
    Z=whos('-file',datafile);
    availableTargetLocs = {Z([Z.bytes]>0 & contains({Z.name},'Target')).name}'
    
    
    
    
    
    
    
    
    %jpsthFigTemplate(datafile,targetLocs,conditions,alignedOn);
    
end