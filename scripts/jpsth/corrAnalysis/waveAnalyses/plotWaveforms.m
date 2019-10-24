
function plotWaveforms(wavformMat,binSize,plotSingleWfFlag,pColor)
    % Convert row = waveform to column = waveform
    wavformMat = wavformMat';
    nPoints = size(wavformMat,1);
    nWaves = size(wavformMat,2);
    % in row form
    x = (1:nPoints).*binSize;
    % plot mean and std
    meanWf = mean(wavformMat,2)';
    stdWf = std(wavformMat,[],2)';
    p = plot(x,meanWf,'LineWidth',2.5,'Color',pColor);
    c = p.Color;
    hold on
    h_var = fill([x fliplr(x)], [meanWf+stdWf fliplr(meanWf-stdWf)],c,'FaceAlpha',0.4,'LineStyle','none');
    set(h_var,'HandleVisibility','off');
    if plotSingleWfFlag
        % individual waveforms
        temp = [wavformMat;nan(1,nWaves)];
        temp = temp(:);
        x = repmat([1:nPoints NaN]',1,nWaves);
        x = x(:);
        ph = plot(x,temp,'Color', c,'LineWidth',0.1);
        ph.Color(4) = 0.1;
    end
    hold off
end