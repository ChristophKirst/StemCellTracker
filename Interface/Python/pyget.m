function v = pyget(varname)
%
% v = pyget(varname)
%
% description:
%   reads a variable form the python workspace and retunrs it in v
%
% See also: pyset

v = py('get', varname);

end

