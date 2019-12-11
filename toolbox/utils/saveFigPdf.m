function [f_h] = saveFigPdf(fn)
    f_h = gcf;
    set(f_h,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',screenposition(3:4),...
        'PaperOrientation','landscape');
    fprintf('Saving figure to: %s...',fn);
    print(fn,'-dpdf','-opengl')
    drawnow
    fprintf('Done\n');
end