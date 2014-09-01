function vig = vignettingFromFlatfield(flatfield, varargin)
%
% vig = vignettingFromFlatfield(flat, param)
% vig = vignettingFromFlatfield(flat, background, param)
%
% description:
%   fits a polynomial of degree 6 and its center to a flat field image
%   the result is an image containing the estimated vignetting
%
% input:
%   flatfield    the flatfield image
%   background   (optional) background image
%   param        paramerter struct with entries
%                optimize.center  optimze center of the vignetting (true)
%                optimize.order   maxim order of the radial polynomial to optimize 2,4,6 (6)
%
% 




end
