function [ unitsByRscSignif ] = categorizeUnitsByRscSignif(rscTbl)
%% Jan 06, 2020
% 2.       What types of neurons contribute to significant r_sc at any point during the trial?
%   a.       Summary plot with SDFs of all such neurons in SEF.
%   b.       Summary plot for FEF.
%   c.       Summary plot for SC.
%   d.       Thomas & Chenchal to meet in 069 to discuss individual SDFs after viewing summary plot.
% 3.       What types of neurons do not contribute to significant r_sc at any point during the trial?
%   a.       Summary plot with SDFs of all such neurons in SEF.
%   b.       Summary plot for FEF.
%   c.       Summary plot forSC.
%   d.       Thomas & Chenchal to meet in 069.
%  
% Yes, and let?s also implement the color density SDF raster plot so we can see variation across all neurons in given groups
%  
%  
% Looking forward to moving this Rich data set forward.
%  
% Me too!
%
%
%%
    % 'Error' = may have both 'ErrorChoice' and 'ErrorTiming'
    %outcomes = {'Correct', 'ErrorChoice', 'ErrorTiming', 'Error'};
    outcomes = {'Correct', 'ErrorChoice', 'ErrorTiming'};
    epochs = {'Baseline', 'Visual', 'PostSaccade', 'PostReward'};
    signifs = [0.05, 0.01];

    unitsFilter.PvalColumn = 'pvalRaw_150ms';
    unitsFilter.Dist = 0;
    % Values redifined in loop below
    unitsFilter.Pval = NaN;
    unitsFilter.Epoch = '';
    unitsFilter.Outcome = '';

    rscTbl.pval = rscTbl.(unitsFilter.PvalColumn);
    % retain only these columns
    rscTbl = rscTbl(:,{'X_unitNum','Y_unitNum','X_area','Y_area','XY_Dist','condition','alignedName','pval'});
    % Convert all NaN distances to Inf to apply distance filter
    rscTbl.XY_Dist(isnan(rscTbl.XY_Dist)) = inf;

    %% For significant and nonSignificant Rsc(s): get unique SEF, FEF, SC units by filtering on different outcomes and epochs
    unitsByRscSignif = table();
    for sig = 1:numel(signifs)
        unitsFilter.Pval = signifs(sig);
        % Add signif. col using values in filt.UseColForPval <= filt.Pval
        rscTbl.sig(rscTbl.pval<=unitsFilter.Pval) = 1;
        for oc = 1:numel(outcomes)
            unitsFilter.Outcome = outcomes{oc};
            for ep = 1:numel(epochs)
                unitsFilter.Epoch = epochs{ep};
                % Split Rsc by other filter criteria
                idxfiltered = ismember(rscTbl.alignedName,unitsFilter.Epoch)...
                    & contains(rscTbl.condition,unitsFilter.Outcome)...
                    & rscTbl.XY_Dist > unitsFilter.Dist;
                unitsByRscSignif = [unitsByRscSignif; getUniqueUnitsByArea(rscTbl(idxfiltered,:),unitsFilter)]; %#ok<*AGROW>
            end
        end
    end
end


%% Sub-functions to get unique units from pair Rsc table
function [unitsByArea] = getUniqueUnitsByArea(rscTblAll,filterCriteria)
   unitsByArea = table();
    signifFlags = [1, 0];
   for s = 1:numel(signifFlags)
    signifFlag = signifFlags(s);
    rscTbl = rscTblAll(rscTblAll.sig == signifFlag,:);
    fns = fieldnames(filterCriteria);
    for fn = 1:numel(fns)
        temp.(['filter_' fns{fn}]) = filterCriteria.(fns{fn});
    end
    temp.filter_IsRscSignificant = signifFlag;
    temp.hasUnitsRscSignificant = 1;
    % get indices for areas
    idxFef = ismember(rscTbl.Y_area,'FEF');
    idxSc = ismember(rscTbl.Y_area,'SC');
    idxSef = ismember(rscTbl.Y_area,'SEF');
    temp.SEF = unique([rscTbl.X_unitNum(idxFef);rscTbl.X_unitNum(idxSc)]);
    temp.FEF = unique(rscTbl.Y_unitNum(idxFef));
    temp.SC = unique(rscTbl.Y_unitNum(idxSc));
    temp.sameArea_SEF_X = unique(rscTbl.X_unitNum(idxSef));
    temp.sameArea_SEF_Y = unique(rscTbl.Y_unitNum(idxSef));    
    % turn into a table
    unitsByArea = [unitsByArea;struct2table(temp,'AsArray',true)];
   
   end
   
   % remove all units from non-significant Rsc that are already accounted
   % for in the significant Rsc
   nonSigRow = find(unitsByArea.filter_IsRscSignificant == 0);
   sigRow = find(unitsByArea.filter_IsRscSignificant == 1);
   temp = unitsByArea(nonSigRow,:);
   temp.SEF{1} = setdiff(unitsByArea.SEF{nonSigRow},unitsByArea.SEF{sigRow});
   temp.FEF{1} = setdiff(unitsByArea.FEF{nonSigRow},unitsByArea.FEF{sigRow});
   temp.SC{1} = setdiff(unitsByArea.SC{nonSigRow},unitsByArea.SC{sigRow});
   temp.sameArea_SEF_X{1} = setdiff(unitsByArea.sameArea_SEF_X{nonSigRow},unitsByArea.sameArea_SEF_X{sigRow});
   temp.sameArea_SEF_Y{1} = setdiff(unitsByArea.sameArea_SEF_Y{nonSigRow},unitsByArea.sameArea_SEF_Y{sigRow});
   temp.hasUnitsRscSignificant = 0; % now we removed all units from significant row
   % add this to units table
   unitsByArea(2,:) = temp;
    
end