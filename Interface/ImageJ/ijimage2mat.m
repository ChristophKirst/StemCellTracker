function data = ijimage2mat(ijplus) 
%
% data = ijimage2mat(ijp)  
%
% description:
%    converts ImageJ ImagePlus image to matalb array in h,w,l coordinates
%
% input:
%    ijplus         ImagePlus java class
%
% output:
%    data           image data
%
% See also: imread_bf


if ~isa(ijplus, 'ij.ImagePlus')
   error('ijimage2mat: expect ij.ImagePlusclass as input image, got: %s!', class(ijplus));
end
   
siz = ijplus.getDimensions()'
data = zeros(siz);

imgs = ijplus.getImageStack.getImageArray;

for i = 1:prod(siz(3:end))
   data(:,:,i) = reshape(imgs(i), siz(1:2));
end

%imagej usually stores data as
%h,w,c,l,t
data = imhwlreshape(data, 'hwclt', 'hwlct');

end