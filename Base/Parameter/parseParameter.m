function param = parseParameter(varargin)
%
% param = parseParameter(varargin)
%
% description:
%    parsed an argument list to a nested parameter struct
%
% input:
%    name, value   pairs, names can contain '.' for nesting
%    param         nested param struct (may be empty struct)
%    []            
%
% output:
%    param         nested param struct

param = varargin2parameter(varargin{:});

end