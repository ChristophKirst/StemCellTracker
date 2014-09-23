function img = imuvwpermute(img, in_format, out_format)
%
% img = imuvwpermute(img, in_format, out_format)
%
% description: 
%     takes a cell array and permutes it from in_format to out_format (in_format = 'uvw' by default)
%
% input:
%     img        input image
%     in_format  (optional) format of input image as characters (imformat(stack))
%     out_format (optional) format of output image, special names: 'matlab' = 'yuw'
%
% output:
%     img       reshaped image
%
% note:
%     'x' = 'u', 'y' = reverse 'v', 'z' = reverse 'w'
%
% See also: impqlpermute

if nargin < 2 || isempty(in_format)
   in_format = imformat(img);
   in_format = strrep(in_format, 'p', 'u');
   in_format = strrep(in_format, 'q', 'v');
   in_format = strrep(in_format, 'l', 'w');
end
if nargin < 3 || isempty(out_format)  % autoconvert to matlab img format
   out_format = 'matlab';
end

if strcmp(in_format, 'matlab')
   in_format = 'yuw';
end
if strcmp(out_format, 'matlab')
   out_format = 'yuw';
end

if length(unique(in_format)) ~= length(in_format)
   error('imuvwpermute: labels appear more than one in input format: %s', in_format)
end
if length(unique(out_format)) ~= length(out_format)
   error('imuvwpermute: labels appear more than one in output format: %s', out_format)
end


%in_format = strrep(in_format, 'x', 'u');
%in_format = strrep(in_format, 'z', 'w');

%out_format = strrep(out_format, 'x', 'p');
%out_format = strrep(out_format, 'z', 'w');


%reverse in y if requested
xposin = strfind(in_format, 'x');
if ~isempty(xposin)
   img = flip(img, xposin);
   in_format(xposin) = 'u';
%else
%   pimg = img;
end
xpos = strfind(out_format, 'x');
if ~isempty(xpos)
   out_format(xpos) = 'u';
end


%reverse in y if requested
yposin = strfind(in_format, 'y');
if ~isempty(yposin)
   img = flip(img, yposin);
   in_format(yposin) = 'v';
%else
%   pimg = img;
end
ypos = strfind(out_format, 'y');
if ~isempty(ypos)
   out_format(ypos) = 'v';
end


zposin = strfind(in_format, 'z');
if ~isempty(zposin)
   img = flip(img, zposin);
   in_format(yposin) = 'w';
%else
%   pimg = pimg;
end
zpos = strfind(out_format, 'z');
if ~isempty(zpos)
   out_format(zpos) = 'w';
end
%in_format
%out_format


%shape = ones(length(out_format),1);
for i = 1:length(out_format)
   sf = strfind(in_format, out_format(i));
   if ~isempty(sf)
      %shape(i) = size(img,sf(1));
   else
      in_format(end+1) = out_format(i);  %#ok<AGROW> % add to the end as this would introduce a new dimensions of size 1 in the original array
   end
end
%shape 
%in_format
%out_format



% find dimensions to squeeze 
for i = 1:length(in_format)
   sf = strfind(out_format, in_format(i));
   if isempty(sf) % removing of dimensions only if it is of size 1 (squeeze)
      if size(img,i) == 1
         out_format = [out_format, in_format(i)]; %#ok<AGROW> %% appending it effectively removes dimensions of size 1
      else
         error('imuvwpermute: dimension %s not specified in output format %of but is non-trivial', in_format(i), out_format);
      end
   end
end
% now out_format and in_format have same size
%in_format
%out_format
%ypos

per = ones(1, length(in_format));
for i = 1:length(in_format)
   per(i) = strfind(in_format, out_format(i));
end
%per

% premute and reshape the img coordinates
img = permute(img, per);

%reverse in y if requested
if ~isempty(xpos)
   xpos = strfind(out_format, 'u');
   img = flip(img, xpos);
end

%reverse in y if requested
if ~isempty(ypos)
   ypos = strfind(out_format, 'v');
   img = flip(img, ypos);
end

%reverse in y if requested
if ~isempty(zpos)
   zpos = strfind(out_format, 'w');
   img = flip(img, zpos);
end


end
   