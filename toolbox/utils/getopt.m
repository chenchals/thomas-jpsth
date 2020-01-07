function x = getopt(argv, options)
% GETOPT checks and returns the arguments of functions
%
% getopt(argv, options) checks that the arguments supplied in the argv cell
% array confrom to the options set in the options cell array. The parsed
% arguments are returned in a structure with appropriately named fields.
%
% The `argv' cell array should be an ordered list of the arguments supplied
% to the function in question. This 1xn argv cell array can be supplied
% directly by varargin.
%
% The `options' cell array should contain a list of possible arguments in a
% fashion similar to GNU getopt(3). An equals ('=') sign at the end of the
% option indicates that this option requires an argument. Lack of an
% equals sign indicates that the respective argument is a flag (i.e. does
% not require an argument). Default values can be specified by enclosing
% the option (switch) and default value in a single cell array. Parameters
% which serve as flags should not specify a default value.
%
% 
% EXAMPLE:
%  Consider a function `foo' defined as foo(varargin). This function takes
%  several arguments with required values (arg1, arg2, arg3). In addition,
%  the function can optionally take a flag (flag1).
%
%  At the beginning of foo, getopt would be called using:
%    getopt(varargin, {'arg1=', 'arg2=', 'arg3=', 'flag1'})
%  A call to foo like foo('arg2', 'x') would cause getopt to return a
%  structure with a field `arg2' whose value was the string 'x' (the
%  structure would also contain a field `flag' as false). A call to
%  foo like foo('arg1', 3, 'flag') would return a structure with a field
%  called `arg1' whose value is the integer 3 and a field `flag' whose
%  value was true.
%
%  Alternatively, if the default values for arg1='x', arg2='y', arg3='z',
%  then they could be specified by a call to getopt using:
%    getopt(varargin, {{'arg1=', 'x'}, {'arg2=', 'y'}, {'arg3=', 'z'}, ...
%    'flag1'})
%

if nargin ~=2
    error('getopt requires exactly two arguments');
end

if isempty(options)
    error('The options cell array must be specified');
end

if ~iscell(options)
    error('The options argument must be a cell array');
end

if ~iscell(argv)
    error('The argv argument must be a cell array');
end

% Construct the output arguments
x = [];

% Check that the contents of options is valid
for i = 1:length(options)
    if iscell(options{i})
        if length(options{i}) ~= 2
            error('Default arguments must be specified using cell arrays whose length is exactly 2');
        end
        if ~ischar(options{i}{1})
            error(['Unknown option format at index ' num2str(i)]);
        end
    elseif ischar(options{i})
        % pass
    else
        error(['Unknown option format at index ' num2str(i)]);
    end
end

% OK, our options array is in the correct format. 

% Set any default values in the structure
for i = 1:length(options)
    if iscell(options{i})
        x.(get_option(options{i}{1})) = options{i}{2};
    elseif ~requires_argument(options{i})
        x.(options{i}) = false;
    end
end

% Now let's parse the arguments list (argv)
skip = false;
for i = 1:length(argv)
    if skip
        skip = false;
        continue;
    end
    found = false;
    for j = 1:length(options)
        if strcmp(get_option(options{j}), argv{i})
            found = true;
            if requires_argument(options{j})
                if i == length(argv)
                    error(['Option ' get_option(options{j}) ' requires an argument']);
                end
                x.(get_option(options{j})) = argv{i+1};
                skip = true;
            else
                % This is a flag
                x.(options{j}) = true;
            end
            break;
        end
    end
    if ~found && ~isempty(argv{i})
        error(['Unknown argument "' argv{i} '"']);
    end
end
end % end of function getopt
    
function requires_arg = requires_argument(x)
    % Returns true if the options x requires an argument
    option = x;
    if iscell(x)
        option = x{1};
    end
    if strcmp(option(length(option)), '=')
        requires_arg = true;
    else
        requires_arg = false;
    end
end

function option = get_option(x)
    % Returns the option (removing any '=' signs)
    option = x;
    if iscell(x)
        option = x{1};
    end
    if strcmp(option(length(option)), '=')
        option = option(1:length(option)-1);
    end
end
