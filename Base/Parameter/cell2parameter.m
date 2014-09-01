function param = cell2parameter(incell)
%
% param = cell2parameter(cell)
%
% description
%    returns a struct with parameters given by varargin
%    paramter names 'xxx.yyy' are parsed to sub structs
%
% input:
%   cell        cell of pairs parname, parvale
% 
% ouput:
%  param        parameter struct
%
% See also: parameter2cell, varargin2cell


param = varargin2parameter(incell{:});
         
end


