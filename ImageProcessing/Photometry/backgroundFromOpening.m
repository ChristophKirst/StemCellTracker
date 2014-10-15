function bkg = backgroundFromOpening(img, varargin)
%
% bkg = backgroundFromOpening(img, param)
%
% description:
%    estimate background using morphological opening
%
% input:
%    imgs    cell array of images, or ImageSource.celldata
%    param   (optional) parameter struct with entries
%            .strel        structural element for opening (strel('disk', 50))
%
% output:
%    bkg     estimated background image 
%
% See also: backgroundFromMin

param = parseParameter(varargin);

if isa(img, 'ImageSource')
   img = img.data(param);
end

se  = getParameter(param, 'strel', 50);

if isnumeric(se)
   se = strel('disk', se);
end

bkg = imopen(img, se);

end