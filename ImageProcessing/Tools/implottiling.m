function implottiling(imgs, varargin)
%
% implottiling(imgs, param)
% implottiling(imgs, titles, param)
%
% descritption:
%    distribute imgs over figure and link thier axes
%
% input:
%    imgs    cellarray of imgs
%    titles  (optional)  cells of title strings
%    param   paramete struct with entries:
%            titles     cell struc of titles ({})
%            tiling     array [max p,max q] of how to tile plane ([] = automatic)
%            link       link zoom on images (true)
%            clf        run clf before plotting (true)
%            clim       color scale, []= automatic, [cmin, cmax] or 'mat2gray' to scale images individually ([])
%            format     'pql' or 'matlab' see note ('pql')
%            other options as in imsubplot
%
% note: to be consistent we assume imgs{p,q} as with pql coordinates
%       in particular imgs{p,q} is left of imgs{p+1,q} and
%                     imgs{p,q} is below of imgs{p,q+1}
%       numbering p is form left lower image to right upper image, 
%       sarting from left to right and then down to up, to use matlab ordering use 'format' option
%
% See also: imsubplot, immontage

titles = {};
if nargin > 1
   if iscellstr(varargin{1})
      titles = varargin{1};
      varargin = varargin(2:end);
   end
end

param = parseParameter(varargin{:});

clim = getParameter(param, 'clim', []);

if ~iscell(imgs)
   % check if we have 3d stack
   if isempty(clim)
      clim = [min(imgs(:)), max(imgs(:))];
   end
   switch imfrmtFormat(imgs)
      case 'XYZ'
         imgs = mat2cell(imgs, size(imgs,1), size(imgs,2), ones(size(imgs,3),1));
         imgs = imgs(:);
      case 'XYZC'
         imgs = mat2cell(imgs, size(imgs,1), size(imgs,2), ones(size(imgs,3),1), size(imgs,4));
         imgs = imgs(:);
      case 'XYCZ'
         imgs = mat2cell(imgs, size(imgs,1), size(imgs,2), size(imgs,3), ones(size(imgs,4),1));  
         imgs = imgs(:);
      otherwise
         imgs = {imgs};
         
         if isempty(clim)
            clim = cscale(imgs);
         end
   end
   %size(imgs)  
   %class(imgs)
else
   if isempty(clim)
      clim = cscale(imgs);
   end
end

if strcmp(clim, 'mat2gray')
   clim = [];
end

if getParameter(param, 'clf', true);
   clf;
end

titles = getParameter(param, 'titles', titles);
tiling = getParameter(param, 'tiling', size(imgs));
ntitles = length(titles(:));
%imgs = imgs'; -> see note


frm = getParameter(param, 'format', 'pql');
if strcmp(frm, 'matlab')
   imgs = imgs';
   imgs = flip(imgs,2);
end


k = 1;
for i = 1:length(imgs(:))
   if ~isempty(imgs{i})
      ax(k) = imsubplot(tiling(1), tiling(2),i, param); %#ok<AGROW>
      k = k + 1;
      if ~isempty(clim)
         implot2d(imgs{i},  'color.scale', clim);
      else
         implot2d(imgs{i});
      end
      if i <= ntitles && ~isempty(titles(i))
         title(titles(i));
      end
   end
end   


if getParameter(param, 'link', true)
   linkaxes(ax, 'xy')
end


end


% helper

function clim = cscale(imgs)
   
   img = imgs{1};
   cmax = double(max(img(:)));
   cmin = double(min(img(:)));

   for i = 2:numel(imgs)
      img = imgs{i};
      cmax = max(cmax, double(max(img(:))));
      cmin = min(cmin, double(min(img(:))));
   end
   
   clim = [cmin, cmax];
end

