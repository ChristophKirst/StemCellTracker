function vign = hivignetting(imgs, shifts, varargin)
%
% imgs = hivignetting(imgs, shifts, varargin)
%
% description:
%    optimize photometrics of images by estimating vignetting (r^6 polynomial) and exposure
%    with help of overlapping pixels
%
% input:
%    imgs    images as cell array of filenames or image data
%    shifts  relative shifts of images
%    param   parameter struct with entries
%            project.filename      project filename (tempname)
%            project.cleanup       cleanup temporary project filename (true)
%            image.filename        temporary image file header (tempname)
%            image.cleanup         cleanup temporary image files (true)
%            options.link          options for which parameters to link
%            options.otimize       paramter to optimize 
%            options.autooptimiser options for autooptimiser ('-n')
%
% output:
%    vign    vignetting image (assumes images are same size)
%
% See also: histitch, hialign

global hitools;
if isempty(hitools)
   hiinitialize();
end

nimgs = length(imgs(:));
if nimgs <= 1
   error('hivignetting: no multiple images to estimate vignetting');
end

param = parseParameter(varargin{:});

pread  = getParameter(param, 'project.read', false);
pclean = getParameter(param, 'project.cleanup', false);
if ~pread
   pclean = false;
end

iclean = getParameter(param, 'image.cleanup', false);

% generate the hugin project and image files

ptofile = hiptogenerate(imgs, shifts, param, 'project.cleanup', false, 'image.cleanup', false);


% optimize photometrics

opts = getParameter(param, 'options.autooptimiser', '-n');
cmd = [hitools('autooptimiser'), ' -o ', ptofile, ' ', opts, ' ', ptofile];
ret = system(cmd);
if ret
   error('hialign: error optimising via command: %s', cmd);
end


% read shifts and clean up

if pread
   shifts = hipto2shifts(ptofile);
   shifts = reshape(shifts, size(imgs));
else
   shifts = ptofile;
end

if pclean
   delete(ptofile);
end

if tclean
   fnlist = strsplit(strtrim(fnlist));
   delete(fnlist{:});
end

end



end