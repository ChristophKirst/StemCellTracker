function img = imrecast(img, clso, varargin)
%
% img = imrecast(data, class, varargin)
%
% description:
%   recast data to class cls
%
% inputs:
%   img     image data
%   class   class name to convert the rescaled image
%   param   (optional) parameter struct with entries
%           rescale    type of rescaling: 'none', 'full' = rescale full data range to new class range, 
%                      [min, max] rescale this data range to new class range ([] = 'none')
%
% output:
%   img     image of new class
%
% See also: imrescale, mat2gray

cls = class(img);
if strcmp(cls, clso)
   return
end

param = parseParameter(varargin{:});

rc = getParameter(param, 'rescale', []);
if isempty(rc)
   rc = 'none';
end

if isnumeric(rc)
   if numel(rc) ~= 2
      error('imrecast: rescale must be [min, max]')
   end
   dmin = rc(1);
   dmax = rc(2);
elseif strcmp(rc, 'full')
   dmin = min(img(:));
   dmax = max(img(:));
end

if strcmp(rc, 'none')
   img =  cast(img, clso);
else
   clsrange = imvaluerange(clso);
   rmin = clsrange(1);
   rmax = clsrange(2); 
   
   rmin = cast(rmin, clso);
   rmax = cast(rmax, clso);

   if dmin == dmax
      img = zeros(size(img), clso) + rmin;
      return
   end

   fac = double(rmax - rmin) / double(dmax - dmin);
   img = cast(double(img - dmin) * fac + double(rmin), clso);
end



