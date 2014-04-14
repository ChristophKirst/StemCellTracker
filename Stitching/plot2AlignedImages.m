function plot2AlignedImages(img1, img2, shift)
%
% plot2AlignedImages(img1, img2, shift)
%

dim = length(shift);
shift1 = zeros(1,dim);
shift2 = zeros(1,dim);

for d = 1:length(shift)
   if shift(d) > 0
      shift2(d) = shift(d);
   else
      shift1(d) = -shift(d);
   end
end

siz = max([size(img1) + shift1; size(img2) + shift2]);
imga = zeros(siz);
imgb = zeros(siz);
imga = imreplace(imga, img1, shift1 + 1);
imgb = imreplace(imgb, img2, shift2 + 1);

imga = imgray2color(imga, 'green');
imgb = imgray2color(imgb, 'magenta');

implot(imga + imgb)

end