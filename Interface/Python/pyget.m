function v = pyget(varname)
%
% v = pyget(varname)
%
% description:
%   reads a variable form the python workspace and retunrs it in v
%
% See also: pyset

v = py('get', varname);
ord = 1:ndims(v);
ord(1) = 2; ord(2) = 1;
v = permute(v, ord);

end

