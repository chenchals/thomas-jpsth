
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
    pltW = 0.26;
    pltH = 0.12;
    gutter = 0.03;
    offsetsX = 0.05:(pltW+gutter*2):1-gutter; % for 3 column-starts
    %offsetsY = [ 0.90 fliplr((0:4).*(pltH+gutter))+0.05]; % for 6 row-starts  
    offsetsY = 0.97-pltH:-(pltH+gutter):gutter; % for 6 row-starts  
    offsetsY(4:end) = offsetsY(4:end) - 0.03;
    pltCount = 0;
    for cols = 1:3
        pos(1) = offsetsX(cols);
        pos(3:4) = [pltW pltH];
        for ros = 1:6
            pos(2) = offsetsY(ros);
            pltCount = pltCount + 1;
            H_plots(pltCount) = axes('parent',H_Figure,'position',pos,'box','on', 'layer','top','Tag',sprintf('H_plot%d',pltCount)); %#ok<AGROW>
        end
    end