function [ binfo , varargout ] = utilIsolateMonkeyBehavior( monkey , binfo , varargin )
%utilIsolateMonkeyBehavior Summary of this function goes here
%   Detailed explanation goes here

if (nargin > 2)
  moves = varargin{1};
else
  moves = zeros(1,length(binfo));
end

if (nargin > 3)
  movesPP = varargin{2};
else
  movesPP = zeros(1,length(binfo));
end

MIN_NUM_TRIALS = 500;

%initializations
sessRemove = [];

%first remove trials from Da/Eu based on presence of SEF data
if ismember({'D'}, monkey)
  sessDa = find(ismember({binfo.monkey}, {'D'}));
  sessRemove = [sessRemove, sessDa(1)];
end
if ismember({'E'}, monkey)
  sessEu = find(ismember({binfo.monkey}, {'E'}));
  sessRemove = [sessRemove, sessEu(1)];
end

binfo(sessRemove) = [];
moves(sessRemove) = [];
movesPP(sessRemove) = [];

%second remove trials based on trial count (MIN_NUM_TRIALS)
idxMonkey = ismember({binfo.monkey}, monkey);
idxNumTrials = ([binfo.num_trials] > MIN_NUM_TRIALS);

binfo = binfo(idxMonkey & idxNumTrials);
moves = moves(idxMonkey & idxNumTrials);
movesPP = movesPP(idxMonkey & idxNumTrials);

if (nargout > 1)
  varargout{1} = moves;
  if (nargout > 2)
    varargout{2} = movesPP;
  end
end

end % util : utilIsolateMonkeyBehavior()

