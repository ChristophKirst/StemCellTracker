function val = getData(obj, varargin)
%
% val = getValue(obj, varargin)
%

val = obj;
var = {varargin{:}}; %#ok<CCAT1>

for i=1:length(var)
   val = [val.(var{i})];
end

end