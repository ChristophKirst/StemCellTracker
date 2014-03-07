function [colmap, coldata] = imlabelcolormap(nlabel, param)
%
% colmap = imlabelcolormap(param)
% colmap = imlabelcolormap(label, param)
%
% description: 
%     generates color map fromo the parameter struct as passed to several plotting routines
%
% input:
%     nlabel    labeled image
%     param     (optional) parameter struct with entries
%               .color.data   color according to this data ([] = random)
%               .color.map    use this colormap (colormap)
%               .color.scale  scale color data (true)
%
% output:
%     colmap    the color map such that cmap(l,:) = colormap(colordata(l))
%     coldata   (optional) color data
%
% See also: colormap, imcolorize, imsurfaceplot

if nargin < 2
   param = [];
end

if isstruct(nlabel)
   param = nlabel;
   nlabel = legnth(param.color.data);
end

cdata  = getParameter(param, {'color', 'data'}, []);
cmap   = getParameter(param, {'color', 'map'}, colormap);
cscale = getParameter(param, {'color', 'scale'}, true);

if isempty(cdata)
   cdata = 1:nlabel;
elseif nlabel ~= length(cdata)
   error('imlabelcolormap: inconsistent color data size!')
end

if ~isempty(cscale) && cscale
   minc = min(cdata);
   maxc = max(cdata);
   if minc ~= maxc
      cdata = (cdata - minc)/(maxc-minc);
   else
      cdata = ones(1, length(cdata));
   end
end

ncm = length(cmap);
colmap = cmap(min(round(cdata * (ncm-1)+1),ncm), :);

if nargout > 1
   coldata = cdata;
end

end
   
