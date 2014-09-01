function pair = alignImagePair(pair, varargin)
%
% pair = alignImagePair(pair, varargin)
%
% description:
%     aligns two images using the AlignmentPair pair and parameters param
%
% input:
%     pair         AlignmentPair class with full images in .from and .to     
%     param        parameter struct with entries
%                  as in align2Images, align2ImagesOnGrid
%
% output:
%     pair         the updated AlignmentPair
%
% See also: align2Images, align2ImagesOnGrid

o = pair.orientation;

if o == 0
   img1 = pair.from;
   img2 = pair.to;
   
   [shift, quality] = align2Images(img1, img2, varargin{:}); 
else
   c = pair.toCell();
   [shift, quality] = align2ImagesOnGrid(c, varargin{:});
end

pair.shift = shift;
pair.quality = quality;

end
