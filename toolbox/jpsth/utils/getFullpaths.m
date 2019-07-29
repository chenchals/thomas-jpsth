function [out] = getFullpaths(pathPattern)
%GETFULLPATHS Get fullpath filenames in the pattern
    d = dir(pathPattern);
    out = strcat({d.folder}', filesep, {d.name}');
end
