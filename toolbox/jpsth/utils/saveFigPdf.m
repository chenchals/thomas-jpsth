function saveFigPdf(fn)
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',screenposition(3:4),...
        'PaperOrientation','landscape');
    fprintf('Saving figure to: %s\n',fn);
    print(fn,'-dpdf','-opengl')
    drawnow
end