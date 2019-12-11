areaPairs ={
    'SEF','SEF'
    'SEF','FEF'
    'SEF','SC'
    'FEF','FEF'
    'FEF','SC'
    'SC','SC'
    'SEF','NSEFN'
    'FEF','NSEFN'
    'NSEFN','NSEFN'
    'SC','NSEFN'
    };
tic
for aa = 1:size(areaPairs,1)
    createSatJpsthDataset(areaPairs{aa,1},  areaPairs{aa,2});  
    toc
end
toc   
    