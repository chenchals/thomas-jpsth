function [ ] = plot_EndptSS_Distr_2(behavInfo, secondSacc )
%plot_EndptSS_Distr Summary of this function goes here
%   varargin : {'D','E'} for monkey

%% filter criteria
monks = {'D','E'};
MIN_NUM_TRIALS = 500;
fx_octant2angle = @(vect) (pi/4)*(double(vect)-1);

%% convert to table for filtering
%note necessary - already (re-)formatted as a table
% behavInfo = struct2table(behavInfo);
% secondSacc = struct2table(secondSacc);

%% Filter sessions by monks
idx = cellfun(@(m) find(strcmp(behavInfo.monkey,m)),monks,'UniformOutput',false);
% remove 1st session for each monk
idx = cellfun(@(x) x(2:end),idx,'UniformOutput',false);
idx = vertcat(idx{:});
behavInfo = behavInfo(idx,:);
secondSacc = secondSacc(idx,:);

%% Filter sessions by min. num trials
idx = find(behavInfo.num_trials>MIN_NUM_TRIALS);
% back to struct for rest of code
behavInfo = table2struct(behavInfo(idx,:));
secondSacc = table2struct(secondSacc(idx,:));

%%
NUM_SESS = size(behavInfo,1);

TGT_ECCEN = 8; %use a consistent eccentricity for plotting

xFinPP = [];
yFinPP = [];
conds = [3, 1];
condColors = {[0 .7 0] [1 0 0]};
condAlphas = [0.1 0.05];
legTxt = {'3 - Fast','1 - Accurate'};
figure();
for c = 1:2
    axesHandle = subplot(1,2,c);
    polaraxes('Units',axesHandle.Units,'Position',axesHandle.Position);
    delete(axesHandle);
    condition = conds(c);
    condColor = condColors{c};
    condAlpha = condAlphas(c);    
    for kk = 1:NUM_SESS        
        %use a consistent target eccentricity
        if (behavInfo(kk).tgt_eccen(100) ~= TGT_ECCEN); continue; end
        
        %index by saccade clipping
        idxClipped = (secondSacc(kk).clipped);
        %index by condition
        %idxCond = (binfo(kk).condition == 3 | binfo(kk).condition == 1);
        %index by condition
        idxCond = (behavInfo(kk).condition == condition);
        
        %index by trial outcome
        idxErr = (behavInfo(kk).err_dir & ~behavInfo(kk).err_time);
        %skip trials with no recorded post-primary saccade
        idxNoPP = (secondSacc(kk).resptime == 0);
        
        %isolate saccade endpoint data
        xfinPP_ = secondSacc(kk).x_fin(idxCond & idxErr & ~idxNoPP & ~idxClipped);
        yfinPP_ = secondSacc(kk).y_fin(idxCond & idxErr & ~idxNoPP & ~idxClipped);
        
        %determine location of singleton relative to absolute right
        th_tgt = fx_octant2angle(behavInfo(kk).tgt_octant((idxCond & idxErr & ~idxNoPP & ~idxClipped)));
        
        %rotate post-primary saccade trajectory according to singleton loc.
        xtmp = cos(2*pi-th_tgt) .* xfinPP_ - sin(2*pi-th_tgt) .* yfinPP_;
        ytmp = sin(2*pi-th_tgt) .* xfinPP_ + cos(2*pi-th_tgt) .* yfinPP_;
        
        xFinPP = cat(2, xFinPP, xtmp);
        yFinPP = cat(2, yFinPP, ytmp);
        
    end%for:session(kk)
    
    
    %% Plotting
    
    %polar distribution of endpoints
    TH_PPSACC = atan2(yFinPP, xFinPP);
    R_PPSACC = sqrt(xFinPP.*xFinPP + yFinPP.*yFinPP);
    
    %polarscatter(TH_PPSACC, R_PPSACC, 40, [.3 .3 .3], 'filled', 'MarkerFaceAlpha',0.3)
    polarscatter(TH_PPSACC, R_PPSACC, 40, condColor, 'filled', 'MarkerFaceAlpha',condAlpha)
    rlim([0 8]); rticks(8); thetaticks([])
    %ppretty([5,5])
    hold on
%     legend(legTxt{c},'Location','northeastoutside')
    
end
%legend(legTxt,'Location','northeastoutside')

end%fxn:plot_EndptSS_Distr()

