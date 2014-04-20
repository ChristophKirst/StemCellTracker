function image = imcolorize(imglab, param)
%
% image = imcolorize(image, param)
%
% description:
%    colorizes the labeled image imglab
%
% input:
%    image      labels, any dimension
%    param      (optional) parameter struct with entries
%               .color.map        use this colormap (colormap)
%               .color.data       color imglab using this value ([] = random)
%               .color.scale      scale color data (true)
%               .color.shuffle    shuffle colors, either 'none', 'shuffle' or 'random' ('shuffle' if color.data = [], 'none' otherwise)
%               .color.background background color
%               
% output:
%    image      colorized labeled rgb image, has one dimension of size 3 more
%
% See also:
%    label2rgb

if nargin < 2
   param = [];
end

cmap  = imlabelcolormap(max(imglab(:)), param);

shuffle = getParameter(param, {'color', 'shuffle'}, []);
if isempty(shuffle)
   shuffle = 'shuffle';
end

back = getParameter(param, {'color', 'background'}, []);
if isempty(back)
    opt = {'k'};
else
    opt = {back};
end

if ~isempty(shuffle)
   opt = [opt, {shuffle}];
end

if ndims(imglab) == 2 %#ok<ISMAT>
   image = label2rgb(imglab, cmap, opt{:});
else
   image = label2rgb3d(imglab, cmap, opt{:});
end

image = double(image) / double(max(image(:)));

end

   

