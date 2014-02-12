function d = imdim(image, dim)
%
% d = imdim(image, dim)
%
% descritpion>
%     returns position of indicated dimension in image
%
% input: 
%     image    the image
%     dim      the dimension 'h', 'w', 'l' = 'z', 'c', 't'
%
% output:
%     d        position of the dimension
%
% See also: imformat

if ~ischar(dim) || isempty(strfind('hwlzct', dim))
   error('imdim: dim must be: h, w, l, z, c, t')
end

switch dim
   case 'z'
      dim = 'l';
   case 'x'    % have to check this for consistency
      dim = 'h';
   case 'y'
      dim = 'w';
end


format = imformat(image);
d = find(format == dim, 1, 'first');

if isempty(d)
   error('imdim: %s not in image with format %s', dim, format)
end

end
