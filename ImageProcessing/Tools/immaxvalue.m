function m = immaxvalue(in)
%
% m = immaxvalue(cls)
%
% description:
%    return max value of an image in, if in is a class name max value of an image of that class
%
% input:
%    in    image data or image class name

m = imvaluerange(in);
m = m(2);

end