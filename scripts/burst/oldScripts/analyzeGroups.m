
    cellInfoTable = T;
    %field mapping:
    session = 'SessionNo';
    sessionValue = 14;
    unitIds = 'SessionCellId';
    channelNo = 'depth';
    
    sessionInfo = cellInfoTable(cellInfoTable.(session)==sessionValue,:);
    
    
    
    selGroup = groups{g};
    %units2Use = eval(selGroup);
    % use these Ids...
    switch selGroup
        case 'units'
            units2Use = arrayfun(@(x) {x},sessionInfo.(unitIds));
            titles = cellfun(@(x) num2str(x,'Unit Id #%02d'),units2Use,'UniformOutput',false);
            analysisType = 'Single Unit';
            pltRows = 5;
            pltCols = 6;
        case 'channelUnits'
            channels  = unique(sessionInfo.(channelNo));
            units2Use = arrayfun(@(x) find(sessionInfo.(channelNo)==x),...
                channels,'UniformOutput',false);
            titles = arrayfun(@(x) num2str(x,'Channel #%02d'),channels,'UniformOutput',false);
            analysisType = 'Pooled units per Channel';
            pltRows = 3;
            pltCols = 6;
        case 'layerUnits'
            units2Use = {
                % upper layer
                cell2mat(arrayfun(@(x) find(sessionInfo.channelNo==x),...
                unique(sessionInfo.channelNo(sessionInfo.upperLayer)),'UniformOutput',false))';
                % lower layer
                cell2mat(arrayfun(@(x) find(sessionInfo.channelNo==x),...
                unique(sessionInfo.channelNo(sessionInfo.lowerLayer)),'UniformOutput',false))';
                };
            titles = {'Upper layers' 'Lower layers'};
            analysisType = 'Pooled units per Layer';
            
            pltRows = 2;
            pltCols = 1;
            
    end






   allBursts = res3.allBursts;
   timeWin = [-500 1500];
    % further burst analysis
    fx_t = @(fn,cellBursts) cellfun(@(x) x.(fn),cellBursts,'UniformOutput',false);     
    for c = 1:size(allBursts,2)
        currBursts = allBursts(:,c);
        allBurstHists(:,c) = BurstUtils.psbh(fx_t('bobT',currBursts),fx_t('eobT',currBursts),timeWin);
    end
   
   
    

    