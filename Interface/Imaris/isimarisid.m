function [b, imaris] = isimarisid(id)
%
% imaris = isimarisid(id)
%
% description:
%    checks if id is a valid Imaris applicarion id
%
% See also: isimaris, imarisinstance

imaris = imarisinstance(id);
b = ~isempty(imaris);

end

