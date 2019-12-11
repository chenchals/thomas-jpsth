
filtCritCellArr = unitSdfsTbl.filterCriteria;

fns2Use = {
    'useEpoch'
    'useOutcome'
    'useRhoPercentile'
    'usePvalForUnit'
    'useMinTrialCount'
    'useRhoUnsigned'
    'rhoUnsignedThresh'
    'rhoPositiveThresh'
    'rhoNegativeThresh'
    'usePvalForPairedUnits'
    };

fltStr = [];
for ii = 1:numel(filtCritCellArr)
    fc = filtCritCellArr{ii};
    tempStr = cell(numel(fns2Use),1);
    for f = 1:numel(fns2Use)
       fn = fns2Use{f};
       fv = fc.(fn);
       if isnumeric(fv)
           tempStr{f} = [fn ' = ' num2str(fv,'%0.3f')];
       elseif ischar(fv)
           tempStr{f} = [fn ' = ' fv];          
       elseif islogical(fv)
           tempStr{f} = [fn ' = ' num2str(fv,'%d')];
           
       end
    end
    fltStr{ii,1} = char(join(tempStr,'; '));
    
end

unique(fltStr,'rows')


