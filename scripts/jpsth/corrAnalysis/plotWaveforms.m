
function plotWaveforms(wavformMat)
  % Each column should be a waveform
  pColor = 'm';
  wavformMat = wavformMat';
  nPoints = size(wavformMat,1);
  nWaves = size(wavformMat,2);
  % in row form
  x = 1:nPoints;
  % plot mean and std
  meanWf = mean(wavformMat,2)';
  stdWf = std(wavformMat,[],2)';
  p = plot(x,meanWf,'LineWidth',1.5,'Color',pColor);
  c = p.Color;
  hold on
  fill([x fliplr(x)], [meanWf+stdWf fliplr(meanWf-stdWf)],c,'FaceAlpha',0.5,'LineStyle','none')
  % individual waveforms
  temp = [wavformMat;nan(1,nWaves)];
  temp = temp(:);
  x = repmat([1:nPoints NaN]',1,nWaves);
  x = x(:);
  ph = plot(x,temp);
  ph.Color = [c(1:3) 0.3];
  hold off
end