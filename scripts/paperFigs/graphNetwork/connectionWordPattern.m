% build a bit pattern of connectivity for 
% functional SEF units
% 
% 1. Filter RscTbl for significant connections to SEF units to FEF and SC
%           keep track of positive and negative Rsc
% 2. Get a unique list of significant SEF units and label their VisMovType
% 3. Parse each connection and build a bit pattern for connection
% 4. Bit Pattern for connection: 
% ________________________________________________________________________
%   Bit name   | bit-value  |  Bit meaning
%   FEFconnect | 2^6-1 = 63 | FEF signif. conn.
%   SCconnect  | 2^5-1 = 31 | SC signif. conn.
%   FEFplus1   | 2^4-1 = 15 | FEF Rsc positive, signif. conn. 1 or more
%   FEFminus1  | 2^3-1 = 7  | FEF Rsc negative, signif. conn. 1 or more
%   SCplus1    | 2^2-1 = 3  | SC Rsc positive, signif. conn. 1 or more
%   SCminus1   | 2^1-1 = 1  | SC Rsc negative, signif. conn. 1 or more
% ------------------------------------------------------------------------
rscTbl = getPairRscTableForHeb();

% remove all rows where rhoRaw_150ms is NaN
rscTbl(isnan(rscTbl.rhoRaw_150ms),:) = [];


% filter rsc table for for SEF to FEF/SC connections only
rscTbl(~ismember(rscTbl.X_area,{'SEF'}),:) = [];
rscTbl(ismember(rscTbl.Y_area,{'SEF'}),:) = [];

% code functional types for X
rscTbl.X_visMovType = repmat({'Other'},size(rscTbl,1),1);
idx = find(abs(rscTbl.X_visGrade) >=2);
rscTbl.X_visMovType(idx) = repmat({'Vis'},numel(idx),1);
idx = find(abs(rscTbl.X_moveGrade) >=2);
rscTbl.X_visMovType(idx) = repmat({'Mov'},numel(idx),1);
idx = find(abs(rscTbl.X_visGrade) >=2 & abs(rscTbl.X_moveGrade) >=2);
rscTbl.X_visMovType(idx) = repmat({'VisMov'},numel(idx),1);
% code: rsc sign
rscTbl.rscSign = sign(rscTbl.rhoRaw_150ms);
% code rsc significance : 
% +1 = rsc positive and significant, -1 = rsc negative and significant
rscTbl.signifRsc = zeros(size(rscTbl,1),1);
idx = find(rscTbl.rhoRaw_150ms > 0 & rscTbl.pvalRaw_150ms <= 0.01);
rscTbl.signifRsc(idx) = ones(numel(idx,1));
idx = find(rscTbl.rhoRaw_150ms < 0 & rscTbl.pvalRaw_150ms <= 0.01);
rscTbl.signifRsc(idx) = ones(numel(idx,1));
% sessions, sef,fef,sc units
sessAll = unique(rscTbl.sess);
sefAll = unique(rscTbl.X_unitNum);
fefAll = unique(rscTbl.Y_unitNum(ismember(rscTbl.Y_area,'FEF')));
scAll = unique(rscTbl.Y_unitNum(ismember(rscTbl.Y_area,'SC')));


%% build the word pattern table for significant connections
sigRscTbl = rscTbl(rscTbl.signifRsc == 1,:);
sefSig = unique(sigRscTbl.X_unitNum);
fefSig = unique(sigRscTbl.Y_unitNum(ismember(sigRscTbl.Y_area,'FEF')));
scSig = unique(sigRscTbl.Y_unitNum(ismember(sigRscTbl.Y_area,'SC')));

% filter table for outcomeEpoch 
epochs = {'PostSaccade'};
outcomes = {'Correct','ErrorChoice','ErrorTiming'};
satConds = {'Fast','Accurate','Fast_Accurate'};

outConnWordPattern = table();
for ep = 1:numel(epochs)
    epoch = epochs{ep};
    for oc = 1:numel(outcomes)
        outcome = outcomes{oc};
        rscs = sigRscTbl(ismember(sigRscTbl.epoch,epoch) & ismember(sigRscTbl.outcome,outcome),:);
        for sc = 1:numel(satConds)
            satCond = satConds{sc};
            if strcmp(satCond,'Fast_Accurate')
                dat = rscs;
            else
                dat = rscs(ismember(rscs.satCondition,satCond),:);
            end
            sefUnitNums = unique(dat.X_unitNum);
            for u = 1:numel(sefUnitNums)
                unitNum = sefUnitNums(u);
                unitPairs = dat(dat.X_unitNum == unitNum,:);
                plusIdx = unitPairs.rscSign == 1;
                minusIdx = unitPairs.rscSign == -1;
                % Count FEF plus/minus pairs
                fefNconnPlus = sum(ismember(unitPairs.Y_area,'FEF') & plusIdx);
                fefNconnMinus = sum(ismember(unitPairs.Y_area,'FEF') & minusIdx);
                % Count SC plus/minus pairs
                scNconnPlus = sum(ismember(unitPairs.Y_area,'SC') & plusIdx);
                scNconnMinus = sum(ismember(unitPairs.Y_area,'SC') & minusIdx);
                % build a table row of output
                wTbl = table();
                wTbl.outcome{1} = outcome;
                wTbl.epoch{1} = epoch;
                wTbl.satCondition{1} = satCond;
                wTbl.sefVisMovType = unique(unitPairs.X_visMovType);
                wTbl.sefUnitNum = unitNum;
                wTbl.nSignifPairs = size(unitPairs,1);
                wTbl.nSignifPlusRsc = sum(plusIdx);
                wTbl.nSignifMinusRsc = sum(minusIdx);
                % parse the pairs for this unit into a bit pattern word
                % initialize bit pattern
                wTbl.FEFconnected = 0; 
                wTbl.SCconnected = 0; 
                wTbl.FEFplus1 = 0;
                wTbl.FEFminus1 = 0;
                wTbl.SCplus1 = 0;
                wTbl.SCminus1 = 0;
                % for FEF pairs
                if fefNconnPlus > 0
                    wTbl.FEFconnected = 1;
                    wTbl.FEFplus1 = 1;
                end
                if fefNconnMinus > 0
                    wTbl.FEFconnected = 1;
                    wTbl.FEFminus1 = 1;
                end
                % for SC pairs
                if scNconnPlus > 0
                    wTbl.SCconnected = 1;
                    wTbl.SCplus1 = 1;
                end
                if scNconnMinus > 0
                    wTbl.SCconnected = 1;
                    wTbl.SCminus1 = 1;
                end
                outConnWordPattern = [outConnWordPattern; wTbl]; %#ok<*AGROW>
                % done for this unit
            end
        end
    end    
end

% add column for 6 bit word
nRows = size(outConnWordPattern,1);
colsFor6bitWord = {'FEFconnected','SCconnected','FEFplus1','FEFminus1','SCplus1','SCminus1'};
colsFor4bitWord = {'FEFplus1','FEFminus1','SCplus1','SCminus1'};
outConnWordPattern.connWord6bit = arrayfun(@(x) sprintf('%d',outConnWordPattern{x,colsFor6bitWord}),(1:nRows)','UniformOutput',false);
outConnWordPattern.connWord4bit = arrayfun(@(x) sprintf('%d',outConnWordPattern{x,colsFor4bitWord}),(1:nRows)','UniformOutput',false);
% convert bit pattern to decimal
outConnWordPattern.connWord6bitDec = arrayfun(@(x) bin2dec(x{1}),outConnWordPattern.connWord6bit);
outConnWordPattern.connWord4bitDec = arrayfun(@(x) bin2dec(x{1}),outConnWordPattern.connWord4bit);

connWordDecToBinPattern = unique(outConnWordPattern(:,{'connWord6bitDec', 'connWord4bitDec','connWord6bit','connWord4bit'}));

% Count different patterns by SEF functional types
epochs = {'PostSaccade'};
outcomes = {'Correct','ErrorChoice','ErrorTiming'};
satConds = {'Fast','Accurate','Fast_Accurate'};
connPatternTbl = outConnWordPattern;
useVars = {'satCondition','sefVisMovType','connWord6bitDec','connWord4bitDec'};
grpVars = useVars; %same order
funcWordWattern = struct();
for oc = 1:numel(outcomes)
    outcome = outcomes{oc};
    idxOutcome = ismember(connPatternTbl.outcome,outcome);
    for ep = 1:numel(epochs)
        epoch = epochs{ep};
        idxEpoch = ismember(connPatternTbl.epoch,epoch);
        for sc = 1:numel(satConds)
            satCond = satConds{sc};
            fn = [outcome epoch satCond];
            idx =  ismember(connPatternTbl.satCondition,satCond) & idxOutcome & idxEpoch;
            temp = connPatternTbl(idx,:);
            res = grpstats(temp(:,useVars),grpVars,{'counts'});
            funcWordWattern.(fn) = res;
        end        
    end
end


              

