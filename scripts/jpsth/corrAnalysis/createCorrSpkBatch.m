areaPairs ={
    'SEF','SEF'
    'SEF','FEF'
    'SEF','SC'
    'FEF','FEF'
    'FEF','SC'
    'SC','SC'
%     'SEF','NSEFN'
%     'FEF','NSEFN'
%     'NSEFN','NSEFN'
%     'SC','NSEFN'
    };
tic
for ii = 1:size(areaPairs,1)
    createCorrSpkDataset(areaPairs{ii,1},areaPairs{ii,2});
    toc
end
toc