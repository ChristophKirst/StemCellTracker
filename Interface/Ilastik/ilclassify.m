function imgseg = ilclassify(ilc, img)
%
% imgseg = ilclassify(img)
%
% description:
%    classifies an imgage img according to the loaded classifier
%
% input:
%    clfile   classifier file as generated with ilastic software
%
% See also: illoadclassifier

if ischar(ilc)
   ilc = illoadclassifier(ilc);
end

if ischar(img)
   if ~isfile(img)
      error('ilclassify: file does not exists!')
   end  
   [~, bi, ei] = fileparts(img);
   img = fullfile(absolutepath(img), [bi, ei]); 
   
   ilc.load_image(img)
else
   imgpy = numpyFromMat(img);
   ilc.image = imgpy;
end

res = ilc.run();
imgseg = squeeze(numpyToMat(res));

end







