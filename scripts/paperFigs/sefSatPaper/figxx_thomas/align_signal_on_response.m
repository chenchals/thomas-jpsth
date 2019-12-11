function [ sig_shift ] = align_signal_on_response( sig , rt )
%align_signal_on_response Summary of this function goes here
%   Detailed explanation goes here

[num_trials, num_samples] = size(sig);

if (num_trials ~= length(rt))
  error('Number of trials in inputs "sig" and "rt" do not match')
end

%make sure the RT data is a float
if ~isfloat(rt)
  rt = double(rt);
end

%initialize re-aligned signal
sig_shift = NaN*ones(num_trials,num_samples);

%loop over all trials

for jj = 1:num_trials
  
  %if response time is zero or NaN, just continue
  if (rt(jj) == 0) || isnan(rt(jj)); continue; end
  
  %shift signal over by response time
  sig_shift(jj,:) = circshift(sig(jj,:), [0,-rt(jj)]);
  
  %remove signal values that were circularly shifted
  if (rt(jj) > 0)
    sig_shift(jj, num_samples-rt(jj)+1:num_samples) = NaN;
  else
    sig_shift(jj, 1:abs(rt(jj))) = NaN;
  end
  
end

end

