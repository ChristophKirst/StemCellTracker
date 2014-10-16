function label = syntheticLabeledImage(isize, nlabel, objsize)
%
% label = syntheticLabeledImage(size, nlabel. objsize)
%
% description:
%     generates a test labeled image
%

label = zeros(isize);
dim = length(isize);

if nargin < 3
   objsize = 10 * ones(1, dim);
else
   objsize = padright(objsize, dim, 'circular');
end

if dim == 2
   for i = 1:nlabel
      wh = round(objsize .* rand(1,2) + 2);
      xy = round((isize-1) .* rand(1,2) + 1);
      label = label + imreplace(zeros(isize), fspecial2('disk', wh, 0) > 0, xy, 'chop', true);
   end
else
   for i = 1:nlabel
      wh = round(objsize .* rand(1,3) + 2);
      xy = round((isize-1) .* rand(1,3) + 1);
      label = label + imreplace(zeros(isize), fspecial3('disk', wh, 0) > 0, xy, 'chop', true);
   end
end

label = bwlabeln(label > 0);

end

      