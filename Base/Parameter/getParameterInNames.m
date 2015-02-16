function [par, pos] = getParameterInNames(param, field, set, varargin)
%
% par = getParameterInNames(param, field, default, set)
%
% input:
%   param     struc of parameters
%   field     the sub field where parameter is stored
%   set       the possible set of names the parameter can assume
%
% output:
%   par       parameter value
%   pos       pos       position in set, 0 if not in set
%
% See also: setParameter, parseParameter


if nargin < 3
   set = {};
end

if isempty(set)
   warning('getParameterInNames: name set empty!');
   par = [];
   pos = 0;
   return
end

if nargin > 3
   default = varargin{1};
else
   default = set{1};
end

par = getParameter(param, field, default);

if isempty(par)
   par = set{1};
   pos = 1;
   return
end

if ischar(par)
   pos = find(ismember(set, par), 1, 'first');
  
   if isempty(pos)
      error('getParameterInNames: unknown name: %s not in %s!', par, var2char(set));
   end
   
elseif isnumeric(par)
   if par > length(set);
      warning('getParameterInNames: parameter id = %d > set length %d, default to: %d!', par, length(set), length(set));
      par = set{end};
      pos = length(set);
   elseif par < 1
      warning('getParameterInNames: parameter id = %d < 1, default to: 1 = %s!', par, set{1});
      par = set{1};
      pos = 1;
   else
      pos = par;
      par = set{pos};
   end
   
else
   error('getParameterInNames: parameter %s not id or name!', var2char(par));   
end

end
   


