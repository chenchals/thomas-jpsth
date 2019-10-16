% load extracted data
dat = load('dataProcessed/analysis/spkCorr/rscSignalNoiseStatic.mat');
fns = fieldnames(dat);
% area pairs to plot (in-order)
pairAreas = {
    'SEF_SEF'
    'SEF_FEF'
    'SEF_SC'
    'SEF_NSEFN'
    'FEF_FEF'
    'FEF_SC'
    'FEF_NSEFN'
    'SC_SC'
    'SC_NSEFN'
    'NSEFN_NSEFN' 
    };
%scatterDataFns = {'
for aa = 1:numel(pairAreas)
     pairArea = pairAreas{aa};
     %scaterData = dat.(pairArea)(:,{
    
    
end


