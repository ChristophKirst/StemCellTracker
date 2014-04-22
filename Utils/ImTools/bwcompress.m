function mm = bwcompress(img, sizei)
%
% mm = bwcompress(img, sizei)
%
% description:
%     use run length encoding on a 2d or 3d bw image input either as array (with
%     values 0, !=0), or as a list of nonzero pixel values (in which case sizei
%     must be supplied).
%
% input:
%   img     bw image array, sizei is irrelevant
%           OR a list of pixels indices of 1's
%   sizei   the size of array as row vector, needed if img input as pixel list.
%
% output:
%   mm      [sizei, (beg end) pix value pairs for intervals where img > 1, col-wise]
% 
% note:
%     For 10x images, this routine uses about 30 uint32 ints per nucleus. There
%     is a matlab bwpack that just bit maps image to unit32 and thus gets a
%     compression of 32.
%
% See also: bwuncompress, bwpack

if min(size(img)) > 1  % assuming input real image
    pixels = find(img);
    sizei = size(img);
else
    if nargin ==1
        error('bwcompress: bad input: pixel list, without sizei\n');
    else
        pixels = img;
    end
end

dim = length(sizei);

npixels = length(pixels);
pixels = sort(pixels);

mm = uint32(zeros(1,2*npixels+2) );  % worst case all regions 1 pixel
mm(1:dim+2) = [dim, sizei, pixels(1)];
last = pixels(1);
ptr = dim+3;
for i = 2:npixels
    if pixels(i) > last +  1;
        %mm = [mm, last, pixels(i)];  %% this took 1000x longer with no
        %prealloc of space.
        mm(ptr) = last;
        mm(ptr+1) = pixels(i);
        ptr = ptr + 2;
        last = pixels(i);
    else
        last = last + 1;
    end
end

mm = [mm(1:(ptr-1)), pixels(end)];

return