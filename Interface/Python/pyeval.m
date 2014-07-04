function varargout = pyeval(cmd)
%
% pyeval(cmd)
% v = pyeval(cmd)
%
% description:
%   evaluates the python code cmd
%
% See also: pyset, pyget

if nargout > 0
   varargout{1} = py('eval', cmd);
else
   py('eval', cmd);
end

end

