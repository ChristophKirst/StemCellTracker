function model = implot3d(image, param)
%
% model = implot3d(image, param)
%
% description:
%    implot3d visualizes 3d volumetirc grayscale or color data in pixel coordinates 
%    using semi transparent slicing
%
% input:
%    image     3d grayscale or color image
%    param     (optional) parameter struct with entries
%              .color.scale    color data scale [cmin, cmax] ([] = automatic = [min(cdata), max(cdata)])
%              .color.alpha    alpha data (image intensity)
%              .renderslice    set the slices to render ('z')
%              .range.p        1x2 p-axis bounds. ( [1 size(data,1)] )
%              .range.q        1x2 q-axis bounds. ( [1 size(data,2)] )
%              .range.l        1x2 l-axis bounds. ( [1 size(data,3)] )
%              .boxratios      1x3 ratio array to scale the axes individually ([1 1 1])
%              .parent         parent axes (gca)
%
% output:
%    model     refernece to the 3d model
%
% note:
%    use interp3 on input data to increase/decrease resolution of data
%
% See also alphamap, colormap, opengl, isosurface

% based on code by Joe Conti, 2004
% modified by C. Kirst 2014 to fit into imtools
%                           in particular assuming pixel coordinates for the data

if nargin < 2
   param = [];
end

if isstruct(image)
    model = image;
else
    model = localGetDefaultModel; 
end

model.cdata = image;

model.parent = getParameter(param, {'parent'}, model.parent);
model.alpha = getParameter(param, {'color', 'alpha'}, model.alpha);
model.cscale = getParameter(param, {'color', 'scale'}, model.cscale);
model.texture = getParameter(param, {'texture'}, model.texture);
model.xdata = getParameter(param, {'range', 'x'}, model.xdata);
model.xdata = getParameter(param, {'range', 'p'}, model.xdata);
model.ydata = getParameter(param, {'range', 'y'}, model.ydata);
model.ydata = getParameter(param, {'range', 'q'}, model.ydata);
model.zdata = getParameter(param, {'range', 'y'}, model.zdata);
model.zdata = getParameter(param, {'range', 'l'}, model.zdata);

rr  = getParameter(param, {'range'}, []);
if ~isempty(rr)
   model.xdata = rr(1:2);         
   model.ydata = rr(3:4);         
   model.zdata = rr(5:6);
end

% define [x,y,z]data
siz = size(model.cdata);
if isempty(model.xdata)
    model.xdata = [0 siz(1)] + 0.5;  % we use pixel coordinates
end
if isempty(model.ydata)
    model.ydata = [0 siz(2)] + 0.5;
end
if isempty(model.zdata)
    model.zdata = [0 siz(3)] + 0.5;
end


boxratios = getParameter(param, {'BoxRatios'}, [1 1 1]);
if ~isempty(boxratios) 
   model.xdata = [model.xdata(1), boxratios(1) * (model.xdata(2)-model.xdata(1)) + model.xdata(1)];
   model.ydata = [model.ydata(1), boxratios(2) * (model.ydata(2)-model.ydata(1)) + model.ydata(1)];
   model.zdata = [model.zdata(1), boxratios(3) * (model.zdata(2)-model.zdata(1)) + model.zdata(1)];
end

if isempty(model.parent)
    model.parent = gca;
end

[model] = local_draw(model);

%set(gca,'YDir','Reverse');
xlabel('p'); ylabel('q'); zlabel('l')

end

%------------------------------------------%
function [model] = localGetDefaultModel

model.cdata = [];
model.alpha = [];
model.xdata = [];
model.ydata = [];
model.zdata = [];
model.parent = [];
model.handles = [];
model.texture = 'xyz';
model.cscale = [];
tag = tempname;
model.tag = ['vol3d_' tag(end-11:end)];

end

%------------------------------------------%
function [model,ax] = local_draw(model)

cdata = model.cdata;
if ndims(cdata) == 3
   cdata = permute(cdata, [2 1 3]);
else
   cdata = permute(cdata, [2 1 3 4]);
end
csiz = size(cdata);

try
   delete(model.handles);
catch
end

ax = model.parent;
%cam_dir = camtarget(ax) - campos(ax);
%[~,ind] = max(abs(cam_dir));

view = zeros(3,1);
switch model.texture
   case 'automatic'
      %view(ind) = 1;
      view(3) = 1; % data is from z stacks
   case 'all'
      view = ones(3,1);
   otherwise
      if ~isempty(strfind(model.texture, 'x'))
         view(1) = 1;
      end
      if ~isempty(strfind(model.texture, 'y'))
         view(2) = 1;
      end
      if ~isempty(strfind(model.texture, 'z'))
         view(3) = 1;
      end
end


opts = {'Parent',ax,'cdatamapping',[],'alphadatamapping',[],'facecolor','texturemap','edgealpha',0,'facealpha','texturemap','tag',model.tag};

if ndims(cdata) > 3
    opts{4} = 'direct';
else
    if isempty(model.cscale)
       cdata = cdata - min(cdata(:));
       cdata = cdata / max(cdata(:));
       caxis(ax, [min(cdata(:)), max(cdata(:))]);
    else
       cdata = cdata - model.cscale(1);
       cdata(cdata < 0) = 0;
       cdata = cdata / (model.cscale(2) - model.cscale(1));
       cdata(cdata> 1) = 1;
       caxis(ax, model.cscale)
    end
    opts{4} = 'scale';
end

if isempty(model.alpha) 
    alpha = cdata;
    if ndims(alpha) > 3
        alpha = sqrt(sum(double(alpha).^2, 4));
        alpha = alpha - min(alpha(:));
        alpha = alpha / max(alpha(:));
    end
    opts{6} = 'none';
else
    alpha = model.alpha;
    if ~isequal(csiz(1:3), size(alpha))
        error('imshow3d: Incorrect size of alpha!');
    end
    opts{6} = 'none';
end

h = findobj(ax,'type','surface','tag',model.tag);
for n = 1:length(h)
  try
     delete(h(n));
  catch
  end
end

handle_ind = 1;

% Create z-slice
if (view(3))    
  x = [model.xdata(1), model.xdata(2); model.xdata(1), model.xdata(2)];
  y = [model.ydata(1), model.ydata(1); model.ydata(2), model.ydata(2)];
  z = [model.zdata(1), model.zdata(1); model.zdata(1), model.zdata(1)];
  diff = model.zdata(2)-model.zdata(1);
  delta = diff/csiz(3);
  for n = 1:csiz(3)

   cslice = squeeze(cdata(:,:,n,:));
   aslice = double(squeeze(alpha(:,:,n)));
   h(handle_ind) = surface(x,y,z,cslice,'alphadata',aslice,opts{:});
   z = z + delta;
   handle_ind = handle_ind + 1;
  end

end


% Create x-slice
if (view(1)) 
  x = [model.xdata(1), model.xdata(1); model.xdata(1), model.xdata(1)];
  y = [model.ydata(1), model.ydata(1); model.ydata(2), model.ydata(2)];
  z = [model.zdata(1), model.zdata(2); model.zdata(1), model.zdata(2)];
  diff = model.xdata(2)-model.xdata(1);
  delta = diff/csiz(2);
  for n = 1:csiz(2)

   cslice = squeeze(cdata(:,n,:,:));
   aslice = double(squeeze(alpha(:,n,:)));
   h(handle_ind) = surface(x,y,z,cslice,'alphadata',aslice,opts{:});
   x = x + delta;
   handle_ind = handle_ind + 1;
  end
end

% Create y-slice
if (view(2))
  x = [model.xdata(1), model.xdata(1); model.xdata(2), model.xdata(2)];
  y = [model.ydata(1), model.ydata(1); model.ydata(1), model.ydata(1)];
  z = [model.zdata(1), model.zdata(2); model.zdata(1), model.zdata(2)];
  diff = model.ydata(2)-model.ydata(1);
  delta = diff/csiz(1);
  for n = 1:csiz(1)

   cslice = squeeze(cdata(n,:,:,:));
   aslice = double(squeeze(alpha(n,:,:)));
   h(handle_ind) = surface(x,y,z,cslice,'alphadata',aslice,opts{:});
   y = y + delta;
   handle_ind = handle_ind + 1;
  end
end

model.handles = h;

end
