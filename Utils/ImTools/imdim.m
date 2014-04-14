function d = imdim(image, dim)
%
% d = imdim(image, dim)
%
% descritpion>
%     returns number of non-color dimensions or if dim is specified the position of indicated dimension in image
%
% input: 
%     image    the image
%     dim      the dimension 'p'='x', 'q'='y', 'l' = 'z', 'c', 't'
%
% output:
%     d        position of the dimension
%
% See also: imformat

if nargin == 1
   format = imformat(image);
   d = length(format);
   d = d - length(strfind(format, 'c'));
else

   if ~ischar(dim) || isempty(strfind('pxqylzct', dim))
      error('imdim: dim must be: p, q, l, x, y, z, c, t')
   end
   
   switch dim
      case 'z'
         dim = 'l';
      case 'x'
         dim = 'p';
      case 'y'
         dim = 'q';
   end
   
   
   format = imformat(image);
   d = find(format == dim, 1, 'first');
   
   if isempty(d)
      error('imdim: %s not in image with format %s', dim, format)
   end

end

end
