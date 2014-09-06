function [cmap, coldata] = imlabelcolormap(nlabel, param)
%
% cmap = imlabelcolormap(param)
% cmap = imlabelcolormap(nlabel, param)
%
% description: 
%     generates a color map from the parameter struct
%
% input:
%     nlabel    (optional) number of lables in labeled image
%     param     (optional) parameter struct with entries
%               .color.map     use this colormap (colormap)
%               .color.data    color label k according using scalar data(k) or color data data(k,:) ([])
%               .color.scale   scale color data (true)
%               .color.shuffle shuffle the colors: [] = 'none', 'shuffle' fixed permutation, 'random' random permuation ('shuffle')
%
% output:
%     colmap    the color map such that cmap(lab,:) = colormap(colordata(lab))
%     coldata   (optional) color data
%
% See also: colormap, imcolormap, imcolorize, imsurfaceplot3d

if nargin < 2
   param = [];
end

if isstruct(nlabel)
   param = nlabel;
end

cdata  = getParameter(param, {'color', 'data'}, []);
cmap = colormap;
cmap = getParameter(param, {'color', 'map'}, cmap);
cscale = getParameter(param, {'color', 'scale'}, true);
cshuffle = getParameter(param, {'color', 'shuffle'}, 'shuffle');
if ~ischar(cshuffle) && cshuffle
   cshuffle = 'shuffle';
else
   cshuffle = 'none';
end
 
if isstruct(nlabel)
   if ~isempty(cdata)
      nlabel = length(cdata);
   else
      nlabel = size(cmap,1);
   end
end


if isempty(cdata)
   cdata = 1:nlabel;
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
else
   cdata(cdata >1) = 1;
   cdata(cdata<0) = 0;
end

if ischar(cmap)
   cmap
   cmap = imcolormap(cmap, nlabel);
end


% shuffle the data if requested
ncm = size(cmap,1);
switch cshuffle
   case 'random'     
      cmap = cmap(randperm(nlabel), :);
   case 'shuffle'
      savedState  = RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', 1));
      cmap = cmap(randperm(nlsa), :);
      RandStream.setGlobalStream(savedState); 
end

cmap = cmap(min(round(cdata * (ncm-1)+1),ncm), :);

if nargout > 1
   coldata = cdata;
end

end
   
