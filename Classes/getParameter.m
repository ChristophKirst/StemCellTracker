function par = getParameter(param, field, varargin)
%
% par = getParameter(param, field, default)
% par = getParameter(param, field1, field2, ..., default)
%
% input:
%   param     struc of parameters
%   field*    the sub field where parameter is stored
%
% output:
%   par       parameter value
%
% See also: setParameter, parseParameter

% initialize 

if nargin > 3
   field = {field varargin{1:end-1}};
   default = varargin{end};
end

if nargin > 2
   default = varargin{1};
end

if nargin <= 2
   default = [];
end

if ~isstruct(param)
   par = default;
   return
end

if ischar(field)
   field = strsplit(field, '.');
end
if ~iscellstr(field)
   par = default;
   return
end

% get parameter

for i = 1:length(field)-1
   if isfield(param, field{i})
      param = param.(field{i});
   else
      par = default;
      return
   end
end

if isfield(param, field{end})
   par = param.(field{end});
else
   par = default;
end



