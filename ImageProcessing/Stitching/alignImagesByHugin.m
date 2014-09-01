function shift = alignImagesByHugin(imgs, varargin)
%
% sh = alignImagesByHugin(imgs, varargin)
%
% description:
%    use the hialign function to aling any images imgs via Hugin interface
%
% input:
%     imgs     images to be alligned as cell array, fist image is anchor
%     param    (optional)  struct with entries as in hialign function
%
% output:
%     shift    the shifts between images (first is origin) in pixel coordinates and pixel units
%
% See also:  hialign

shift = hialign(imgs, varargin{:});

end

