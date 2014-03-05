function inv = iminvert(image)
%
% inv = iminvert(image)
%
% description:
%     inverts a grayscale image
%
% input:
%     image    grayscale image to invert
%
% output:
%     inv      inverted image

inv = max(image(:)) - image;

end