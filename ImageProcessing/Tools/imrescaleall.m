function imgs = imrescaleall(imgs, varargin)
%
% imgs = imrescaleall(imgs, varargin)
%
% description:
%   rescale all images in cell array imgs commonly 
%
% inputs:
%   data    cell of image data
%   param   (optional) parameter struct with entries
%           class       class to convert the rescaled images to ([] = class(data))
%           data.min    minimal value for overall data ([] = min(all data(:) in cell))
%           data.max    maximal value for overall data ([] = max(all data(:) in cell))
%           rescale.max max value of rescaled data ([] = max value for class, or 1 for double/single)
%           rescale.min min value of rescaled data ([] = min value for class, or 0 for double/single)
%
% See also: mat2gray

if ~iscell(imgs)
   error('imrescaleall: expects cell of images as input!');
end

if isempty(imgs)
   return
end

param = parseParameter(varargin{:});

dmin = getParameter(param, 'data.min', []);
dmax = getParameter(param, 'data.max', []);

clso = getParameter(param, 'class', class(imgs{1}));
clsrange = imvaluerange(clso);

rmin = getParameter(param, 'rescale.min', clsrange(1));
rmax = getParameter(param, 'rescale.max', clsrange(2)); 

rmin = cast(rmin, clso);
rmax = cast(rmax, clso);

nimgs = numel(imgs);


if isempty(dmin)
   dmin = Inf;
   for i = 1:nimgs
      data = imgs{i};
      dmin = min(dmin, min(data(:)));
   end
end

if isempty(dmax)
   dmax = -Inf;
   for i = 1:nimgs
      data = imgs{i};
      dmax = max(dmax, max(data(:)));
   end
end

if dmin == dmax
   for i = 1:nimgs
      imgs{i} = zeros(size(imgs{i}), clso) + rmin;
   end
   return
end

fac = double(rmax - rmin) / double(dmax - dmin);
for i = 1:nimgs
   imgs{i} = cast(double(imgs{i} - dmin) * fac + double(rmin), clso);
end

end



