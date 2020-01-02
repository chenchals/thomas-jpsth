function [f_h] = saveFigPdf(fn)
%SAVEFIGPDF saves current figure in landscape mode. 
%   To be proper *must* use getFigureHandle to create the base figure
% see also NEWFIGURE
%
    orient = 'landscape';
    f_h = gcf;
    set(f_h,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',screenposition(3:4),...
        'PaperOrientation',orient);
    fprintf('Saving figure to: %s...',fn);
    print(fn,'-dpdf','-opengl')
    drawnow
    fprintf('Done\n');
end
