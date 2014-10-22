function pyimport(varargin)
%
% pyimport(varargin)
%
% description:
%    imports a set of variables to the matlab workspace from python
%

   for i=1:numel(varargin)
		tmp = py('get', varargin{i});
      assignin('base', varargin{i}, tmp);
   end
end
