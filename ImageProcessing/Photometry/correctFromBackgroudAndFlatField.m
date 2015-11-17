function img = correctFromBackgroudAndFlatField(img, bkg, flt, varargin)
%
% img = correctFromBackgroudAndFlatField(img, bkg, flt)
%
% description:
%    corrects image via (img - bkg) / (flt - bkg)
%
% input:
%    imgs    cell array of images, or ImageSource.celldata
%    param   (optional) parameter struct with entries
%            intensity.min     minimal intensity after background substraction
%            division.offset   add this to the denominator for division
%
% output:
%    bkg     estimated background image 
%
% note: 
%    Reference: Fundamentals of Light Microscopy and Electronic Imaging, p 421
% 
% See also: backgroundFromMin

param = parseParameter(varargin);

imin = getParameter(param, 'intensity.min', 10* eps);
doff = getParameter(param, 'division.offset', 0);

noCell = false;
if ~iscell(img)
   img = {img};
   noCell = true;
end

div =  (flt - bkg) + doff;
div(div < imin) = imin;
%min(div(:))

for i = 1:length(img)
   ic = (img{i} - bkg);
   ic(ic < imin) = imin;
   %min(ic(:))
   img{i} = ic ./ div;
end

if noCell
   img = img{1};
end

end