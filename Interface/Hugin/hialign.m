function shifts = hialign(imgs, varargin)
%
% shifts = hialign(imgs, param)
% ptofile = hialign(imgs, param)
%
% description:
%    align images imgs using cpfind and autooptimiser from Hugin
%    cpfind finds control/key points, 
%    autooptimiser aligns images
%
% input:
%    imgs    images as cell array of filenames or image data
%    param   parameter struct with entries
%            project.filename      project filename (tempname)
%            project.cleanup       cleanup temporary project filename (true)
%            project.read          parse pto file and read shifts (true)
%            image.filename        temporary image filename header (tempname)
%            image.cleanup         cleanup temporary image files (true)
%            image.unlink          unlink images for individual lenses (true)
%            options.ptogen        options for pto_gen ('')
%            options.ptovar        options for pto_var ('--opt TrX,TrY')
%            options.cpfind        options for cpfind ('--fullscale --celeste')
%            options.autooptimiser options for autooptimiser ('-n')
%
% output:
%    shifts  relative shifts of images in pixel coordinates and pixel units
%    ptofile for project.read == false the pto filename
%
% See also: histitch

global hitools;
if isempty(hitools)
   hiinitialize();
end

nimgs = length(imgs(:));
if nimgs == 0
   error('hialign: no images to align');
end

param = parseParameter(varargin{:});

pread  = getParameter(param, 'project.read', true);
pclean = getParameter(param, 'project.cleanup', true);
if ~pread
   pclean = false;
end

iclean = getParameter(param, 'image.cleanup', true);

% generate the hugin project and image files

[ptofile, ifilelist] = hiptogenerate(imgs, 'options.ptovar', '--opt TrX,TrY ', param, 'project.read', false, 'image.cleanup', false);

% find control points

%opts = getParameter(param, 'options.cpfind', '--multirow --celeste');
%opts = getParameter(param, 'options.cpfind', '--fullscale --celeste');
%opts = getParameter(param, 'options.cpfind', '--fullscale');
opts = getParameter(param, 'options.cpfind', '--fullscale --sieve1size 1000 --sieve1height 100 --sieve1width 100 --minmatches 2 --sieve2size 5');

cmd = [hitools('cpfind'), ' -o ', ptofile, ' ', opts, ' ', ptofile];
ret = system(cmd);
if ret
   error('hialign: error finding key points via command: %s', cmd);
end


% optimize positions

opts = getParameter(param, 'options.autooptimiser', '-n');
cmd = [hitools('autooptimiser'), ' -o ', ptofile, ' ', opts, ' ', ptofile];
ret = system(cmd);
if ret
   error('hialign: error optimization via command: %s', cmd);
end


% read shifts and clean up

if pread || iclean
   pto = hiparsepto(ptofile);
else
   pto = ptofile;
end

if pread
   shifts = hipto2shifts(pto);
   shifts = reshape(shifts, size(imgs));
else
   shifts = ptofile;
end

%ls
%ifilelist
%iclean

if pclean
   delete(ptofile);
end

if iclean
   delete(ifilelist{:});
end

end
