function rimage = imhwlreshape(image, in_format, out_format)
%
% image = imhwlreshape(image, in_format, out_format)
%
% description: 
%     takes a stack of gray images and moves them into a hwl format
%
% input:
%     image      h x w x l image stack
%     in_format  (optional) format of ipnut image as characters (imformat(stack))
%     out_format (optional) format of output image  ('hwltc')
%
% output:
%     rimage      reshaped image
%
% See also: imformat

if nargin < 2 || isempty(in_format)
   in_format = imformat(image);
end
if nargin < 3 || isempty(out_format)  % autoconvert to matlab image format
   out_format = 'matlab';
end

if strcmp(in_format, 'matlab')
   in_format = 'hwclt';
end
if strcmp(out_format, 'matlab')
   out_format = 'hwclt';
end
%if length(in_format) ~= ndims(image) % matlab is inconssitent if last dimension is 1 
%   error('imhwlreshape: input image and format inconsistent!')
%end
if length(unique(in_format)) ~= length(in_format)
   error('imhwlreshape: labels appear more than one in input format: %s', in_format)
end
if length(unique(out_format)) ~= length(out_format)
   error('imhwlreshape: labels appear more than one in output format: %s', out_format)
end


%find output shape and possible additional dimensions
shape = ones(length(out_format),1);
k = 1;
for i = 1:length(out_format)
   sf = strfind(in_format, out_format(i));
   if ~isempty(sf)
      shape(i) = size(image,sf(1));
      of(k) = out_format(i);  %#ok<AGROW>
      k = k + 1;
   end
end

% find positions
k = 1;
for i = 1:length(in_format)
   sf = strfind(of, in_format(i));
   if isempty(sf) % removing of dimensions only if it is of size 1 (squeeze)
      if size(image,i) == 1
         siz = size(image); siz(i) = [];
         image = reshape(image, siz);
      else
         error('imhwlreshape: dimension %s are not specified in format: %s', in_format(i), format);
      end
   else
      per(k) = sf(1); %#ok<AGROW>
      k = k + 1;
   end
end

% premute and reshape the image coordinates
rimage = reshape(permute(image, per), shape(:)');

end
   