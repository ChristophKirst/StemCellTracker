function img = histitch(imgs, shifts, varargin)
%
% img = histitch(imgs, shifts, param)
%
% description:
%    stitch images imgs using alignments shifts via the enblend tool
%    enblend stitching is based on stiching the image in deifferent levels of Laplacian pyramid:
%    i.e. blending is done on differnet scales for different spatial frequency components
%
% input:
%    imgs    images as cell array or image file names in which case the images must have an alpha channel representing the mask (see note)
%    shifts  shifts of iondividual cells as returned by alignment routines
%    param   parameter struct with entries
%            output.filename    filename of stitched output image (tempname)
%            output.read        read the final stiched image from temporary file (true)
%            output.cleanup     remove output files from disk (true)
%            temporary.filename filename header for temporary tiled image files (tempname)
%            temporary.cleanup  remove temporary image tile files from disk (true)
%            options            enblend options as text ('')
%
% output:
%    img     stiched image
%
% note: an alpha channel can be used to ignore certain regions in one image in the overlap region
%
% See also: hialign, hiinitialize

global hitools;
if isempty(hitools)
   hiinitialize();
end

param = parseParameter(varargin{:});

fn = getParameter(param, 'output.filename', tempname);
[~, ~, fext]  = fileparts(fn);
if isempty(fext)
   fext = '.tif';
   fn = [fn fext];
end

fnt = getParameter(param, 'temporary.filename', tempname);
[ftpath, ftname, ftext]  = fileparts(fnt);
if isempty(ftext)
   ftext = '.tif';
end

nimgs = length(imgs(:));
if nimgs == 0
   error('histitch: no images to stitch');
end

oread =  getParameter(param, 'output.read', true);
oclean = getParameter(param, 'output.cleanup', true);
if ~oread
   oclean = false;
end
tclean = getParameter(param, 'temporary.cleanup', true);
opts   = getParameter(param, 'options', '');

% check images
if iscellstr(imgs)
   
   % list of image filenames
   fnlist = '';
   for i = 1:nimgs
      if ~isfile(imgs{i})
         error('histitch: image file %s does not exists!', imgs{i});
      else
         fnlist = [fnslist, ' ', imgs{i}];
      end
   end
   
else % list of images

   % check dims
   dim = ndims(imgs{1});
   if dim ~=2 && dim ~=3
      error('histitch: images should be gray or rgb, found %g channels!', dim);
   end
   for i = 1:nimgs
      if ~isequal(dim, ndims(imgs{i}))
         error('histitch: image %g has different dimension %g ~= %g', i, dim, ndims(imgs{i}));
      end
   end
   
   %check for trivial stiching 
   if (nimgs == 1)
      img = imgs{1};
      if ~oread
         imwrite_tiff(img, fn);
      end
      if ~tclean
         imwrite_tiff(img, fullfile(ftpath, [ftname num2str0(1, 4) ftext]));
      end
      return
   end
   
   
   %write images to temporary files and create image filename list

   isize = cellfun(@size, imgs, 'UniformOutput', false);
   ashift = absoluteShiftsAndSize(shifts, isize);

   %resolution = 150;
   resolution = 1;

   
   % if images are double we need to rescale them to 0 - 1
   if isequal(class(imgs{1}), 'double')
      imgs = imrescaleall(imgs);
   end
   
   mxval = immaxvalue(class(imgs{1}));
   
   fnlist = '';
   for i = 1:nimgs
      imgo = imgs{i};
 
      if dim==2
         imgo = cat(3, imgo, mxval * ones(size(imgo), 'like', imgo));
      else
         imgo = cat(3, imgo, mxval * ones(size(imgo), 'like', imgo));
      end
   
      %imgo = imrescale(imgo, 'class', 'uint16');
      fni = fullfile(ftpath, [ftname num2str0(i, 4) ftext]);
      fnlist = [fnlist ' ' fni]; %#ok<AGROW>

      % x,y coordinates (exchange from pq)
      sh = ashift{i};
      imwrite_tiff(imgo, fni, 'XPosition', sh(2) / resolution, 'YPosition', sh(1)/resolution, ...
                              'XResolution', resolution, 'YResolution', resolution, ...
                              'ExtraSamples', Tiff.ExtraSamples.UnassociatedAlpha);
   end
end


% enblend

ret = system([hitools('enblend') ' ', opts , ' -o ', fn, ' ', fnlist]);
if ret
   error('histitch: error using enblend, via command: %s', [hitools('enblend') ' ', opts , ' -o ', fn, ' ', fnlist]);
end


% read and clean up

if oread
   img = imread(fn);
   if size(img,3) == 2
      img = img(:,:,1);
   else
      img = img(:,:,1:3);
   end
else
   img = fn;
end

if oclean
   delete(fn);
end

if tclean
   fnlist = strsplit(strtrim(fnlist));
   delete(fnlist{:});
end

end