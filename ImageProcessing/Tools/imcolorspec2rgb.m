function rgb = imcolorspec2rgb(col)
%
% rgb = imcolorspec2rgb(col)
%
% description: converts col to a rgb format
%
% input:
%    col    rgb color, matlab color name
%
% output:
%    rgb    rgb color
%
% See also: cspecchk


if ischar(col)
   switch col
      case {'y', 'yellow'}
         rgb = [1 1 0];
      case {'m', 'magenta'}
         rgb = [1 0 1];
      case {'c', 'cyan'}
         rgb = [0 1 1];
      case {'r', 'red'}
         rgb = [1 0 0];
      case {'g', 'green'}
         rgb = [0 1 0];
	   case {'b', 'blue'}
         rgb = [0 0 1];
	   case {'w', 'white'}
         rgb = [1 1 1];
      case {'k', 'black'}
         rgb = [0 0 0];
      otherwise
         error('imcolor2rgb: color %s not valid!', col)
   end
else
   rgb = col;
end

end
      