function param = setParameter(varargin)
%
% param = setParameter(varargin)
%
% description
%    returns a struct with parameters given by varargin
%    paramter names 'xxx.yyy' are parased to sub structs
%
% input:
%   varargin    parameter paris parname, parvale
% 
% ouput:
%  param        parameter struct
%
% See also: getParameter, parseParameter, cell2parameter

param = cell2parameter(varargin{:});