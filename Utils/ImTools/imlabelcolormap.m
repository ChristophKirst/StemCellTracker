function [colmap, coldata] = imlabelcolormap(nlabel, param)
%
% colmap = imlabelcolormap(param)
% colmap = imlabelcolormap(label, param)
%
% description: 
%     generates a color map from the parameter struct
%
% input:
%     nlabel    number of lables in labeled image
%     param     (optional) parameter struct with entries
%               .color.data    color according to this data ([] = randperm(nlabel))
%               .color.map     use this colormap (colormap)
%               .color.scale   scale color data (true)
%
% output:
%     colmap    the color map such that cmap(l,:) = colormap(colordata(l))
%     coldata   (optional) color data
%
% See also: colormap, imcolorize, imsurfaceplot3d

if nargin < 2
   param = [];
end

if isstruct(nlabel)
   param = nlabel;
   nlabel = length(param.color.data);
end

cdata  = getParameter(param, {'color', 'data'}, []);
cmap = colormap;
cmap   = getParameter(param, {'color', 'map'}, cmap);
cscale = getParameter(param, {'color', 'scale'}, true);

if isempty(cdata)
   %cdata = 1:nlabel;
   
   % shuffle the data -> use same seed
   savedState  = RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', 1));
   cdata = randperm(nlabel);
   RandStream.setGlobalStream(savedState);
   
elseif nlabel ~= length(cdata)
   error('imlabelcolormap: inconsistent color data size!')
end
cdata = double(cdata);

if ~isempty(cscale) && cscale
   minc = min(cdata);
   maxc = max(cdata);
   if minc ~= maxc
      cdata = (cdata - minc)/(maxc-minc);
   else
      cdata = ones(1, length(cdata));
   end
end

ncm = size(cmap,1);
colmap = cmap(min(round(cdata * (ncm-1)+1),ncm), :);

if nargout > 1
   coldata = cdata;
end

end
   
