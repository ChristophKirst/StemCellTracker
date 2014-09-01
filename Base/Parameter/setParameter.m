function param = setParameter(varargin)
%
% param = setParameter(varargin)
% param = setParameter(param, varargin)
%
% description
%    returns a struct with parameters given by varargin
%    paramter names 'xxx.yyy' are parased to sub structs
%
% input:
%   varargin    parameter pairs parname, parvale
%   param       (optional) parameter struct
% 
% ouput:
%  param        parameter struct
%
% See also: getParameter, parseParameter, cell2parameter

param = varargin2parameter(varargin{:});

end