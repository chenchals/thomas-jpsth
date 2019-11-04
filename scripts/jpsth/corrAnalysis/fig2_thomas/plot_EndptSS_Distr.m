function [ ] = plot_EndptSS_Distr( binfo , movesPP , varargin )
%plot_EndptSS_Distr Summary of this function goes here
%   Detailed explanation goes here

args = getopt(varargin, {{'monkey=',{'D','E'}}});

%[binfo, ~, movesPP] = utilIsolateMonkeyBehavior(binfo, zeros(1,length(binfo)), movesPP, args.monkey);
[ binfo , ~, movesPP ] = utilIsolateMonkeyBehavior( args.monkey , binfo , zeros(1,length(binfo)), movesPP );
NUM_SESS = length(binfo);

TGT_ECCEN = 8; %use a consistent eccentricity for plotting

xFinPP = [];
yFinPP = [];
conds = [3, 1];
condColors = {[0 1 0] [1 0 0]};
condAlphas = [0.2 0.1];
legTxt = {'Condition 3','Condition 1'};
figure();
for c = 1:2
    axesHandle(c) = subplot(1,2,c);
    polarAxesHandle(c) = polaraxes('Units',axesHandle(c).Units,'Position',axesHandle(c).Position);
    delete(axesHandle(c));
    condition = conds(c);
    condColor = condColors{c};
    condAlpha = condAlphas(c);
    for kk = 1:NUM_SESS
        
        %use a consistent target eccentricity
        if (binfo(kk).tgt_eccen(100) ~= TGT_ECCEN); continue; end
        
        %index by saccade clipping
        idxClipped = (movesPP(kk).clipped);
        %index by condition
        %idxCond = (binfo(kk).condition == 3 | binfo(kk).condition == 1);
        %index by condition
        idxCond = (binfo(kk).condition == condition);
        
        %index by trial outcome
        idxErr = (binfo(kk).err_dir & ~binfo(kk).err_time);
        %skip trials with no recorded post-primary saccade
        idxNoPP = (movesPP(kk).resptime == 0);
        
        %isolate saccade endpoint data
        xfinPP_ = movesPP(kk).x_fin(idxCond & idxErr & ~idxNoPP & ~idxClipped);
        yfinPP_ = movesPP(kk).y_fin(idxCond & idxErr & ~idxNoPP & ~idxClipped);
        
        %determine location of singleton relative to absolute right
        th_tgt = convert_tgt_octant_to_angle(binfo(kk).tgt_octant((idxCond & idxErr & ~idxNoPP & ~idxClipped)));
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
    rlim([0 8]); thetaticks([])
    %ppretty([5,5])
    hold on
    legend(legTxt{c},'Location','northeastoutside')

end
%legend(legTxt,'Location','northeastoutside')

end%fxn:plot_EndptSS_Distr()
