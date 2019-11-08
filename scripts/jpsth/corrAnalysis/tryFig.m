
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