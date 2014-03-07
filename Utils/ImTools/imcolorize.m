function image = imcolorize(label, param)
%
% image = imcolorize(image, param)
%
% description:
%    colorizes the label image
%
% input:
%    image      matrix of labels
%    param      (optional) parameter struct with entries
%               .color.data       color according to this data ([] = random)
%               .color.map        use this colormap (colormap)
%               .color.scale      scale color data (true)
%               .color.shuffle    shuffle colors
%               .color.background background color
%               
% output:
%    image      colorized labeled image
%
% See also:
%    label2rgb

if nargin < 2
   param = [];
end

cmap  = imlabelcolormap(max(label(:)), param);

shuffle = getParameter(param,  {'color', 'shuffle'}, []);

if isempty(shuffle) && isempty(getParameter(param,  {'color', 'data'}, []))
   shuffle = 'shuffle';
end

back = getParameter(param,  {'color', 'background'}, []);

if isempty(back)
    opt = {'k'};
else
    opt = {back};
end
if ~isempty(shuffle)
   opt = [opt, {shuffle}];
end

if ndims(label) == 2 %#ok<ISMAT>
   image = label2rgb(label, cmap, opt{:});
else
   image = label2rgb3d(label, cmap, opt{:});
end

image = double(image) / double(max(image(:)));

end

   

