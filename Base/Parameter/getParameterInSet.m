function [par, pos] = getParameterInSet(param, field, set, varargin)
%
% par = getParameterInSet(param, field, default, set)
%
% input:
%   param     struc of parameters
%   field     the sub field where parameter is stored
%   set       the possible set of values the parameter can assume
%
% output:
%   par       parameter value
%   pos       pos       position in set, 0 if not in set
%
% See also: setParameter, parseParameter


if nargin < 3
   set = {};
end

if nargin > 3
   default = varargin{1};
else
   if isempty(set)
      default = [];
   else
      default = set{1};
   end
end

par = getParameter(param, field, default);
pos = find(ismember(set, par));

if pos == 0
   par = default;
end

end
   


