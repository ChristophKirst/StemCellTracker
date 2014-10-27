function varargout = pyprint(var)
%
% pyprint(var)
%
% description:
%   runs print var in python
%
% See also: pyeval

py('eval', ['print ' var]);

end

