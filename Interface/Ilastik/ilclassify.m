function imgseg = ilclassify(img)
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


if ischar(img)
   if ~isfile(img)
      error('ilclassify: file does not exists!')
   end  
   [~, bi, ei] = fileparts(img);
   img = fullfile(absolutepath(img), [bi, ei]); 
   
   py('eval', ['ilc.load_image(' img ')'])
else
   py('set', 'img', img)
   py('eval', 'ilc.image = img')
end


py('eval', 'res = ilc.run()')
imgseg = py('get', 'res');
imgseg = squeeze(imgseg);

%imf = imformat(imgseg)
%imfn = strrep(imf, 'q', 'x'); imfn = strrep(imfn, 'p', 'q')
%imgseg = impqlpermute(imgseg, imf, imfn);

end







