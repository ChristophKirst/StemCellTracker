function img = imrescale(img, varargin)
%
% img = imrescale(data, varargin)
%
% description:
%   rescale data into class cls
%
% inputs:
%   data    image data
%   param   (optional) parameter struct with entries
%           class       class to convert the rescaled image to ([] = class(data))
%           data.min    minimal value for overall data ([] = min(data(:)))
%           data.max    maximal value for overall data ([] = max(data(:)))
%           rescale.max max value of rescaled data ([] = max value for class, or 1 for double/single)
%           rescale.min min value of rescaled data ([] = min value for class, or 0 for double/single)
%
% See also: mat2gray

param = parseParameter(varargin{:});

dmin = getParameter(param, 'data.min', min(img(:)));
dmax = getParameter(param, 'data.max', max(img(:)));

clso = getParameter(param, 'class', class(img));
clsrange = imvaluerange(clso);

rmin = getParameter(param, 'rescale.min', clsrange(1));
rmax = getParameter(param, 'rescale.max', clsrange(2)); 

rmin = cast(rmin, clso);
rmax = cast(rmax, clso);

if dmin == dmax
   img = zeros(size(img), clso) + rmin;
   return
end

fac = double(rmax - rmin) / double(dmax - dmin);
img = cast(double(img - dmin) * fac + double(rmin), clso);

end



