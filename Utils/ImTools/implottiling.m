function implottiling(images, tiling, titles)
%
% imtile(images, tiling, titles)
% imtile(images, titles)
%
% descritption:
%    distribute images over figure and link thier axes
%
% input:
%    images    cellarray of images
%    tiling    (optional)  [m, n]  the tiling in rows and columns
%    titles    (optional)  cells of title strings
%
% See also: immontage, imsubplot

if ~iscell(images)
   % check if we have 3d stack
   switch imformat(images)
      case 'pql'
         images = mat2cell(images, size(images,1), size(images,2), ones(size(images,3),1));
      case 'pqlc'
         images = mat2cell(images, size(images,1), size(images,2), ones(size(images,3),1), size(images,4));
      case 'pqcl'
         images = mat2cell(images, size(images,1), size(images,2), size(images,3), ones(size(images,4),1));  
      otherwise
         images = {images};  
   end
   
   if nargin < 2
      tl = ceil(sqrt(length(images)));
      tl = [tl ceil(length(images)/tl)];
      im = cell(tl);
      im(1:length(images)) = images(:);
      images = im;
   end
   
end

if nargin < 2
   tiling = size(images);
   titles = {};
end
if nargin < 3
   if iscell(tiling)
      titles = tiling;
      tiling = size(images);
   else
      titles = {};
   end
end

ntitles = length(titles(:));
images = images';

k = 1;
for i = 1:length(images(:))
   if ~isempty(images{i})
      ax(k) = imsubplot(tiling(1),tiling(2),i); %#ok<AGROW>
      k = k + 1;
      %size(images{i})
      
      implot(images{i});
      if i <= ntitles && ~isempty(titles(i))
         title(titles(i));
      end
   end
end   
linkaxes(ax, 'xy')

end