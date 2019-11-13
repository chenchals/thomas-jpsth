%% spike distribution plot B
set(0,'units','pixels');
set(0,'defaulttextfontsize',6,...
    'defaulttextfontname','Arial',...
    'defaultaxesfontsize',6,...
    'defaultaxeslinewidth',0.05);
margin = 10; %pixels
%ss=get(0,'ScreenSize');
% optimized for this size on my macbookpro
ss = [1 1 1680 1050];
ss = [20 20 1500 990];
FigPos=[margin margin ss(3)-(2*margin) ss(4)-(2*margin)];
%Main figure window
H_Figure=figure('Position',FigPos,...
    'color',[1 1 1],'numbertitle','off','renderer','painters',...
    'renderermode','manual','menubar','none',...
    'Tag','H_Figure');
orient landscape
set(H_Figure,'units','normalized')
ss = get(0,'ScreenSize');
aspectRatio = ss(3)/ss(4);

 nAreas = 3;

    nRows = nAreas*2;
    if nRows == 6
        pltH = 0.12;
        pltW = 0.26;
        gutter = 0.03;
        offsetsX = 0.05:(pltW+gutter*2):1-gutter; % for 3 column-starts
        offsetsY = 0.95-pltH:-(pltH+gutter):gutter; % for 6 row-starts
        offsetsY(4:end) = offsetsY(4:end) - gutter;
    elseif nRows == 4
        pltH = 0.18;
        pltW = 0.275;
        gutter = 0.04;
        offsetsX = 0.05:(pltW+gutter):1-gutter; % for 3 column-starts
        offsetsY = 0.95-pltH:-(pltH+gutter):gutter; % for 4 row-starts
        offsetsY(3:end) = offsetsY(3:end) - gutter;
    end   
    
    pltCount = 0;
    for cols = 1:3
        pltW2 = (pltW-gutter)/2;
        for ros = 1:nRows
            for subp = 1:2 
                % (1) percentage of pairs vs binned-spikecorr, 
                % (2) spike corr vs binned-distance
                pos(1) = offsetsX(cols) + (pltW2+gutter/1.5)*(subp-1);
                pos(3:4) = [pltW2 pltH];
                pos(2) = offsetsY(ros);
                pltCount = pltCount + 1;
                H_plots(pltCount) = axes('parent',H_Figure,'position',pos,'box','on', 'layer','top','Tag',sprintf('H_plot%d',pltCount)); %#ok<AGROW>
            end
        end
    end

%% spike distribution plot A
ss = [20 20 1500 990];
margin = 10;
FigPos=[margin margin ss(3)-(2*margin) ss(4)-(2*margin)];
%Main figure window
H_Figure=figure('Position',FigPos,...
    'color',[1 1 1],'numbertitle','off','renderer','painters',...
    'renderermode','manual','menubar','none',...
    'Tag','H_Figure');
orient landscape
set(H_Figure,'units','normalized')
ss = get(0,'ScreenSize');
aspectRatio = ss(3)/ss(4);

    nRows = 6;
    if nRows == 6
        pltH = 0.12;
        pltW = 0.26;
        gutter = 0.03;
        offsetsX = 0.05:(pltW+gutter*2):1-gutter; % for 3 column-starts
        offsetsY = 0.95-pltH:-(pltH+gutter):gutter; % for 6 row-starts
        offsetsY(4:end) = offsetsY(4:end) - gutter;
    elseif nRows == 4
        pltH = 0.18;
        pltW = 0.275;
        gutter = 0.04;
        offsetsX = 0.05:(pltW+gutter):1-gutter; % for 3 column-starts
        offsetsY = 0.95-pltH:-(pltH+gutter):gutter; % for 4 row-starts
        offsetsY(3:end) = offsetsY(3:end) - gutter;
    end
    pltCount = 0;
    for cols = 1:3
        pos(1) = offsetsX(cols);
        pos(3:4) = [pltW pltH];
        for ros = 1:nRows
            pos(2) = offsetsY(ros);
            pltCount = pltCount + 1;
            H_plots(pltCount) = axes('parent',H_Figure,'position',pos,'box','on', 'layer','top','Tag',sprintf('H_plot%d',pltCount)); %#ok<AGROW>
        end
    end
    
    %% for spike corr comparision among areas
    % 3 rows : correct, errorChoice, errorTiming
    % 4 cols : baseline, visual, postSaccade, postReward
    ss = [20 20 1500 990];
    margin = 10;
    FigPos=[margin margin ss(3)-(2*margin) ss(4)-(2*margin)];
    %Main figure window
    H_Figure=figure('Position',FigPos,...
        'color',[1 1 1],'numbertitle','off','renderer','painters',...
        'renderermode','manual','menubar','none',...
        'Tag','H_Figure');
    orient landscape
    set(H_Figure,'units','normalized')
    ss = get(0,'ScreenSize');
    aspectRatio = ss(3)/ss(4);
    
    nRows = 3;
    nCols = 4;
    
    pltH = 0.25;
    pltW = 0.20;
    gutter = 0.04;
    offsetsX = 0.06:(pltW+gutter):1-gutter; % for 4 column-starts
    offsetsY = 0.95-pltH:-(pltH+gutter+0.01):gutter; % for 3 row-starts
    
    pltCount = 0;
    for cols = 1:nCols
        pos(1) = offsetsX(cols);
        pos(3:4) = [pltW pltH];
        for ros = 1:nRows
            pos(2) = offsetsY(ros);
            pltCount = pltCount + 1;
            H_plots(pltCount) = axes('parent',H_Figure,'position',pos,'box','on', 'layer','top','Tag',sprintf('H_plot%d',pltCount)); %#ok<AGROW>
        end
    end
    

    