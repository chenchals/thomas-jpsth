function [ idx_poor_isolation ] = identify_trials_poor_isolation_SAT( ninfo , num_trials , varargin )
%identify_trials_poor_isolation_SAT Summary of this function goes here
%   Detailed explanation goes here

args = getopt(varargin, {{'task=','SAT'}});

if strcmp(args.task, 'SAT')
  field_ = 'trRemSAT';
elseif strcmp(args.task, 'MG')
  field_ = 'trRemMG';
else
  error('Input "task" not recognized')
end

if isempty(ninfo.(field_)) %no trials to remove due to poor isolation
  idx_poor_isolation = false(1,num_trials);
elseif (ninfo.(field_)(1) == 9999) %remove ALL trials from consideration
  idx_poor_isolation = true(1,num_trials);
else %sepcified interval of trials to remove
  idx_poor_isolation = false(1,num_trials);
  idx_poor_isolation(ninfo.(field_)(1) : ninfo.(field_)(2)) = true;
end

end%util:identify_trials_poor_isolation_SAT()

