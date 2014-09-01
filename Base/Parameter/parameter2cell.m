function cellpar = parameter2cell(param)
%
% cellpar = parameter2cell(param)
%
% description:
%    converts a nested parameter struct to a cell of the form {'name1', val1, ...}
%
% input:
%    param   paramete struct
%
% output:
%    cellpar cellstr with pairs name value, nested structs are converted to dotted names
%
% See also: cell2parameter, varargin2parameter

if nargin == 0 || isempty(param)
   cellpar = {};
   return
end

cellpar = generateCells('', param);

end

function cpar = generateCells(rootname, param)
   cpar = {};
   fns = fieldnames(param);
   for i = 1:length(fns)
      val = param.(fns{i});
      if ~isempty(rootname) 
         fns{i} = [rootname '.' fns{i}];
      end
      
      if isstruct(val)
         cpar = [cpar, generateCells(fns{i}, val)]; %#ok<AGROW>
      else
         cpar = [cpar, fns(i), {val}]; %#ok<AGROW>
      end
   end
end