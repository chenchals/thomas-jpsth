areaPairs ={
    'SEF','SEF'
    'SEF','FEF'
    'SEF','SC'
    'SEF','NSEFN'
    'FEF','FEF'
    'FEF','SC'
    'FEF','NSEFN'
    'SC','SC'
    'SC','NSEFN'
    'NSEFN','NSEFN'
    };
tic
for aa = 1:size(areaPairs,1)
    createSatJpsthDataset(areaPairs{aa,1},  areaPairs{aa,2});  
    toc
end
toc   
    