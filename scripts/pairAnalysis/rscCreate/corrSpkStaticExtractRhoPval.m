%%oFilename = 'dataProcessed/analysis/11-18-2019/spkCorr/summary/spkCorrAllPairsStaticRhoPval.mat';
spkCorrFile = 'dataProcessed/analysis/11-18-2019/spkCorr/summary/spkCorrAllPairsStaticNew.mat';
areaPairFields = {
    'SEF_SEF' 
    'SEF_FEF'    
    'SEF_SC'     
    'FEF_FEF'    
    'FEF_SC'     
    'SC_SC'      
    };
corrDatFields = {
    'Pair_UID'                 
    'X_unitNum'       
    'Y_unitNum'       
    'X_area'          
    'Y_area'          
    'XY_Dist'         
    'condition'          
    'alignedName'  
    'rhoRaw_50ms'
    'pvalRaw_50ms'
    'rhoRaw_150ms'
    'pvalRaw_150ms'
    'rhoRaw_200ms'
    'pvalRaw_200ms'
    };

spkCorrStatic = struct();
for ii = 1:numel(areaPairFields)
    ap = areaPairFields{ii};
    temp = load(spkCorrFile,ap);
    temp = temp.(ap);
    spkCorrStatic.(ap) = temp(:,corrDatFields);
end


save(oFilename,'-v7.3','-struct','spkCorrStatic');

