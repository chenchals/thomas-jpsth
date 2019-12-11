function [ sdf , varargout ] = compute_spike_density_fxn( spike_times )
%compute_spike_density_fxn Summary of this function goes here
%   Detailed explanation goes here

%make sure our input is of the correct format
if ~iscell(spike_times)
  error('Input "spike_times" should be of class "cell"')
end

%initialize the excitatory post-synaptic potential
tau_d = 20; tau_g = 1;
epsp = @(x) exp(-x/tau_d) .* (1 - exp(-x/tau_g));
epsp_conv = epsp(transpose(linspace(0,199,200)));
epsp_conv = epsp_conv * 1000/sum(epsp_conv);

%times (in msec) of trial start and end
TRIAL_START = -3500;
TRIAL_END = 2500;

%get number of trials and smaples based on input
num_trials = length(spike_times);
num_samples = TRIAL_END - TRIAL_START + 1;

%initialize output spike density function (and spike trains)
sdf = single(NaN*ones(num_trials, num_samples));
train = false(num_trials, num_samples);

%loop over all trials and calculate single-trial functions
for kk = 1:num_trials
  
  %get single-trial spike train
  train(kk,spike_times{kk}) = true;
  
  %compute the convolution
  temp = conv(single(train(kk,:)), epsp_conv, 'full');
  sdf(kk,:) = temp(1:num_samples);
  
end

%if desired, return the spike trains
if (nargout > 0)
  varargout{1} = train;
end

end

