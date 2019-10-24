function [wavWidths] = getWaveformWidths(wavforms)
    % given a cell-array of cell-arrays (matrix of waveforms for each trial) :-), 
    % find width by interpolating with a cubic-spline at 10x sampling
    % example: load waveformfile
    % load('Unit_0001.mat)
    % ww = wf.wavSearch -> is {cell_array_of_wavefroms_for_each_trial}
    
    fx_spkWid = @(spks,magFactor) ...
        arrayfun(@(s) (find(spks(s,:)==max(spks(s,:),[],2),1) ...
                    - find(spks(s,:)==min(spks(s,:),[],2),1))/magFactor,...
                      (1:size(spks,1))');
    wavWidths = cell(numel(wavforms),1);
    interpolateAt = 10; % 10x
    for ro = 1:numel(wavforms)
        w = wavforms{ro};
        if isempty(w)
            continue;
        end      
        x = size(w{1},2);
        x = (1:x)';
        xq = (1:1/interpolateAt:max(x))';
        % there could be trials with no spikes
        notEmptyIdx = find(~cellfun(@isempty,w));
        % interpolated waveforms for trials with spikes
        wq = cellfun(@(s) interp1(x,s',xq,'pchip')',w(notEmptyIdx),'UniformOutput',false);
        % waveform widths in sampling space
        wqWideTemp = cellfun(@(s) fx_spkWid(s,interpolateAt),wq,'UniformOutput',false);
        % put back widths for correct trials
        wqWide = cell(size(w,1),1);
        wqWide(notEmptyIdx) = wqWideTemp;
        wavWidths{ro} = wqWide;
        
    end
end

