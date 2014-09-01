function imgs = hiphotometricoptimization(imgs, shifts, varargin)
%
% imgs = hiphotometricoptimization(imgs, shifts, varargin)
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
%            project.read          parse pto file and read shifts (true)
%            image.filename        temporary image file header (tempname)
%            image.cleanup         cleanup temporary image files (true)


%            options.ptogen        options for pto_gen ('')
%            opttionsptovar        options for pto_var ('--opt TrX,TrY --set v=10')
%            options.cpfind        options for cpfind ('--multirow --celeste')
%            options.autooptimiser options for autooptimiser ('-n')
%
% output:
%    shifts  relative shifts of images in pixel coordinates and pixel units
%    ptofile for project.read == false the pto filename
%
% See also: histitch, hialign

global hitools;
if isempty(hitools)
   hiinitialize();
end

nimgs = length(imgs(:));
if nimgs <= 1
   error('hiphotometricoptimization: no multiple images to photometric optimize');
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