function data = ijimage2mat(ijplus) 
%
% data = ijimage2mat(ijp)  
%
% description:
%    converts ImageJ ImagePlus image to matalb array
%
% input:
%    ijplus         ImagePlus java class
%
% output:
%    data           image data
%
% See also: imreadBF


if ~isa(ijplus, 'ij.ImagePlus')
   error('ijimage2mat: expect ij.ImagePlus class as input image, got: %s!', class(ijplus));
end
   
siz = ijplus.getDimensions()';
data = zeros(siz);

imgs = ijplus.getImageStack.getImageArray;

for i = 1:prod(siz(3:end))
   data(:,:,i) = reshape(imgs(i), siz(1:2));
end

%imagej usually stores data as
%y,x,c,l,t
data = impqlpermute(data, 'yxczt', 'xyzct');

end