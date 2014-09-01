function m = imvaluerange(in)
%
% m = imvaluerange(cls)
%
% description:
%    returns [min, max] value of an image in, if in is a class name [min, max] value of an image of that class
%
% input:
%    in    image data or image class name


if ischar(in)
   cls = in;
else
   cls = class(in);
end


switch cls
   case {'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32', 'uint64'}
      m = [intmin(cls), intmax(cls)];
   case {'single', 'double'}
      m = [0,1];
   otherwise
      error('imvaluerange: class %s no supported image data format!', cls);
end

end