function imgres = correctResponse(img, varargin)
%
% [segments, distances] = segmentByPropagation(label, param)
%
% description: 
%   passes teh image through a response function described by param
%
% input:
%    img             instensity image to be transformed
%    param           parameter struct
%                    .vignetting.cente
%
% output:
%

param = parseParameter(varargin);

si = size(img);

ep = getParameter(param, 'exposure', 1);

vc = getParameter(param, 'vignetting.center', si /2);
vp = getParameter(param, 'vignetting.parameter', [1,0,0,0]);

rt = getParameter(param, 'response.type', 1); 
rp = getParameter(param, 'response.parameter', [0,0,0,0,0,0]);

iv = getParameter(param, 'inverse', 0);

img = mexResponseTransform(img, [], [], [], [], ep, vc, vp, rt, rp, iv);


end

   