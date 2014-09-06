function m = imminvalue(in)
%
% m = imminvalue(cls)
%
% description:
%    return min value of an image in, if in is a class name main value of an image of that class
%
% input:
%    in    image data or image class name

m = imvaluerange(in);
m = m(1);

end