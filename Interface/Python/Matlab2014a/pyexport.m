function pyexport(varargin)
%
% pyexport(varargin)
%
% description:
%      exprots a set of variables from the matlab workspace to python
%
	for i=1:numel(varargin)
		py('set', varargin{i}, evalin('base',varargin{i}));
	end
end
