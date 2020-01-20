
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
filterEpochs = {'PostSaccade'} ;%unique(categorizedUnitsTbl.filter_Epoch); %{'PostSaccade'} ;%
filterOutcomes = {'Correct'}; %unique(categorizedUnitsTbl.filter_Outcome); %{'Correct'}; %
filterPvals = 0.05;%unique(categorizedUnitsTbl.filter_Pval); %0.05;%

%% 
for oc = 1:numel(filterOutcomes)
    useOutcome = filterOutcomes{oc};
    for ep = 1:numel(filterEpochs)
        useEpoch = filterEpochs{ep};
        for pv = 1:numel(filterPvals)
            usePval = filterPvals(pv);
            str = sprintf('%s_%s_%s',useOutcome,useEpoch,num2str(usePval*100,'%03d'));
            [satSdfsTbl.(str),satSdfsImageTbl.(str)] = plotSatSdfsHeatmapByArea(categorizedUnitsTbl,useOutcome,useEpoch,usePval);
        end        
    end
end
