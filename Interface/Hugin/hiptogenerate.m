function [ptofile, ifilelist] = hiptogenerate(imgs, varargin)
%
% ptofile = hiptogenerate(imgs, param)
% ptofile = hiptogenerate(imgs, shifts, param)
% [ptofile, ifilelist] = hiptogenerate(...)
%
% description:
%    generates pto file from images imgs
%
% input:
%    imgs      images as cell array of filenames or image data
%    shifts    (optional) shifts of images in pixel coordinates
%    param     parameter struct with entries
%              project.filename      project filename (tempname)
%              project.read          parse project file (false)
%              project.cleanup       delete project file (false)
%              image.filename        temporary image file name header (tempname)
%              image.cleanup         cleanup temporary image files (false)
%              image.unlink          unlink images for individual lenses (true)
%              options.ptogen        options for pto_gen ('')
%              options.ptovar        additional options for pto_var ('')
%
% output:
%    ptofile   pto filename
%    ifilelist image file list
%
% See also: histitch, hialign

global hitools;
if isempty(hitools)
   hiinitialize();
end

nargin
if nargin > 1 && iscell(varargin{1})
   shifts = varargin{1};
   varargin = varargin(2:end);
else
   shifts = {};
end
%shifts

param = parseParameter(varargin{:});

ptofile = getParameter(param, 'project.filename', tempname);
[ptopath, ptoname, ptoext]  = fileparts(ptofile);
if isempty(ptoext)
   ptoext = '.pto';
end
ptofile = fullfile(ptopath, [ptoname ptoext]);

pread  = getParameter(param, 'project.read', false);
pclean = getParameter(param, 'project.cleanup', false);
if ~pread
   pclean = false;
end


ifile = getParameter(param, 'image.filename', tempname);
[ipath, iname, iext]  = fileparts(ifile);
if isempty(iext)
   iext = '.tif';
end
iclean = getParameter(param, 'image.cleanup', false);
 

nimgs = length(imgs(:));
if nimgs == 0
   error('hiptogenerate: no images for pto file!');
end


% check images
if iscellstr(imgs)
   
   % list of image filenames
   ifilelist = '';
   for i = 1:nimgs
      if ~isfile(imgs{i})
         error('hiptogenerate: image file %s does not exists!', imgs{i});
      else
         ifilelist = [fnslist, ' ', imgs{i}];
      end
   end
   
else % list of images

   % check dims
   dim = ndims(imgs{1});
   if dim ~=2 && dim ~=3
      error('hiptogenerate: images should be gray or rgb, found %g channels!', dim);
   end
   for i = 1:nimgs
      if ~isequal(dim, ndims(imgs{i}))
         error('hiptogenerate: image %g has different dimension %g ~= %g', i, dim, ndims(imgs{i}));
      end
   end
      
   %write images to temporary files and create image filename list
   %resolution = 150;
   resolution = 1;

   ifilelist = '';
   for i = 1:nimgs
      imgo = imgs{i};

      fni = fullfile(ipath, [iname num2str0(i, 4) iext]);
      ifilelist = [ifilelist ' ' fni]; %#ok<AGROW>

      % x,y coordinates (exchange from pq)
      imwrite_tiff(imgo, fni, 'XResolution', resolution, 'YResolution', resolution);
   end
end

% generate pto file

opts = getParameter(param, 'options.ptogen', '');
cmd = [hitools('pto_gen'), ' ', opts , ' -o ', ptofile, ' ', ifilelist];
ret = system(cmd);
if ret
   error('hiptogenerate: error generating pto file via command: %s', cmd);
end

% set lens parameter

% we dont optimize for TrZ=0
% => field of view has to be adjusted for each lens to match image size

pto = hiparsepto(ptofile);

setv = '--set ';
v0 = 1; % degree horizontal field of view
zref = pto(1).w / 2 / tan(v0/2 * 2*pi/360);
for i = 1:nimgs
   v = 2 * atan2(pto(i).w /2, zref) /2/pi*360;
   setv = [setv, 'v', num2str(i-1), '=' num2str(v), ',']; %#ok<AGROW>
end
setv = setv(1:end-1);

cmd = [hitools('pto_var'), ' -o ', ptofile, ' ',  setv, ' ', ptofile];
ret = system(cmd);
if ret
   error('hiptogenerate: error changing pto file via command: %s', cmd);
end


% write shifts
if ~isempty(shifts)
   pto = hiparsepto(ptofile);
   isizes = hipto2sizes(pto);  
   pto = hishifts2pto(shifts, isizes);
   
   setTr = '--set ';
   
   for i = 1:length(pto)
      setTr = [setTr, 'TrX', num2str(i-1), '=', num2str(pto(i).TrX), ',']; %#ok<AGROW>
      setTr = [setTr, 'TrY', num2str(i-1), '=', num2str(pto(i).TrY), ',']; %#ok<AGROW>
      setTr = [setTr, 'TrZ', num2str(i-1), '=', num2str(pto(i).TrZ), ',']; %#ok<AGROW>
   end
   
   cmd = [hitools('pto_var'), ' -o ', ptofile, ' ', setTr, ' ', ptofile];
   ret = system(cmd);
   if ret
      error('hiptogenerate: error changing pto file via command: %s', cmd);
   end
end

% unlink images to get individual lenses for each image / unlink vignetting ?
ulink = getParameter(param, 'image.unlink', true);
if isempty(ulink)
   ulink = false;
elseif iscell(ulink)
   ulinkvars = ulink;
   ulink = true;
else
   ulinkvars = {'v', 'Ra', 'Rb', 'Rc', 'Rd', 'Re', 'a', 'b', 'c', 'd', 'e', 'g', 't', 'Va', 'Vb', 'Vc', 'Vd', 'Vx', 'Vy'};
end
   
if ulink
   ulink = '--unlink ';
   for u = 1:length(ulinkvars)
      for i = 0:nimgs-1
         ulink = [ulink, ulinkvars{u}, num2str(i), ',']; %#ok<AGROW>
      end
   end
   ulink = ulink(1:end-1);
   
   cmd = [hitools('pto_var'), ' -o ', ptofile, ' ', ulink, ' ', ptofile];
   ret = system(cmd);
   if ret
      error('hiptogenerate: error changing pto file via command: %s', cmd);
   end
end

% additional pto_var options
opts = getParameter(param, 'options.ptovar', '');
if ~isempty(opts)
   cmd = [hitools('pto_var'), ' -o ', ptofile, ' ', opts, ' ', ptofile];
   ret = system(cmd);
   if ret
      error('hiptogenerate: error changing pto file via command: %s', cmd);
   end
end

% parse pto file
if pread
   ptofile = hiparsepto(ptofile);
end

% clean up
if pclean
   delete(ptofile);
end

ifilelist = strsplit(strtrim(ifilelist));

if iclean
   delete(ifilelist{:});
end

end
