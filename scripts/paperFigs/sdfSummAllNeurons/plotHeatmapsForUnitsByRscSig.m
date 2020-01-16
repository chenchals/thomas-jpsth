
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
%% Categorize units by filtering on R_sc significance
categorizedUnitsTbl = categorizeUnitsByRscSignif();

% Filter criteria that can be used to categorize units by Rsc
filterEpochs = unique(categorizedUnitsTbl.filter_Epoch);
filterOutcomes = unique(categorizedUnitsTbl.filter_Outcome);
filterPvals = unique(categorizedUnitsTbl.filter_Pval);

%% 
for oc = 1:numel(filterOutcomes)
    useOutcome = filterOutcomes{oc};
    for ep = 1:numel(filterEpochs)
        useEpoch = filterEpochs{ep};
        for pv = 1:numel(filterPvals)
            usePval = filterPvals(pv);
            [sadSdfsTbl,satSdfsImageTbl] = plotSatSdfsHeatmapByArea(categorizedUnitsTbl,useOutcome,useEpoch,usePval);
        end        
    end
end