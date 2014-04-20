function model = implot3d(img, varargin)
%
% model = implot3d(img, varargin)
%
% description:
%    implot3d visualizes 3d volumetirc grayscale or color data in pixel coordinates 
%    using semi transparent slicing
%
% input:
%    img       3d grayscale or color image
%    param     (optional) parameter struct with entries
%              .color.scale    color data scale [cmin, cmax] ([] = automatic = [min(cdata), max(cdata)])
%              .color.alpha    alpha data (img intensity)
%              .renderslice    set the slices to render ('z')
%              .range.p        1x2 p-axis bounds. ( [1 size(data,1)] )
%              .range.q        1x2 q-axis bounds. ( [1 size(data,2)] )
%              .range.l        1x2 l-axis bounds. ( [1 size(data,3)] )
%              .boxratios      1x3 ratio array to scale the axes individually ([1 1 1])
%              .parent         parent axes (gca)
%
% output:
%    model     reference to the 3d model
%
% note:
%    use interp3 on input data to increase/decrease resolution of data
%
% See also alphamap, colormap, opengl, isosurface

% based on code by Joe Conti, 2004
% modified by C. Kirst 2014 to fit into imtools

if nargin < 2
   param = [];
elseif nargin > 3
   param = setParameter(varargin{:});
else
   param = varargin{1};
end

if isstruct(img)
    model = img;
else
    model = localGetDefaultModel; 
end

model.cdata = img;

model.parent = getParameter(param, {'parent'}, model.parent);
model.alpha = getParameter(param, {'color', 'alpha'}, model.alpha);
model.cscale = getParameter(param, {'color', 'scale'}, model.cscale);
model.texture = getParameter(param, {'texture'}, model.texture);
model.xrange = getParameter(param, {'range', 'x'}, model.xrange);
model.xrange = getParameter(param, {'range', 'p'}, model.xrange);
model.yrange = getParameter(param, {'range', 'y'}, model.yrange);
model.yrange = getParameter(param, {'range', 'q'}, model.yrange);
model.zrange = getParameter(param, {'range', 'y'}, model.zrange);
model.zrange = getParameter(param, {'range', 'l'}, model.zrange);

rr  = getParameter(param, {'range'}, []);
if ~isempty(rr)
   model.xrange = rr(1:2);         
   model.yrange = rr(3:4);         
   model.zrange = rr(5:6);
end

% define [x,y,z] data
siz = size(model.cdata);
if isempty(model.xrange)
    model.xrange = [0 siz(1)] + 0.5;
end
if isempty(model.yrange)
    model.yrange = [0 siz(2)] + 0.5;
end
if isempty(model.zrange)
    model.zrange = [0 siz(3)] + 0.5;
end


boxratios = getParameter(param, {'BoxRatios'}, [1 1 1]);
if ~isempty(boxratios) 
   model.xrange = [model.xrange(1), boxratios(1) * (model.xrange(2)-model.xrange(1)) + model.xrange(1)];
   model.yrange = [model.yrange(1), boxratios(2) * (model.yrange(2)-model.yrange(1)) + model.yrange(1)];
   model.zrange = [model.zrange(1), boxratios(3) * (model.zrange(2)-model.zrange(1)) + model.zrange(1)];
end

if isempty(model.parent)
    model.parent = gca;
end

[model] = local_draw(model);

%set(gca,'YDir','Reverse');
xlabel('p'); ylabel('q'); zlabel('l')
xlim(model.xrange); ylim(model.yrange); zlim(model.zrange);

end

%------------------------------------------%
function [model] = localGetDefaultModel

model.cdata = [];
model.alpha = [];
model.xrange = [];
model.yrange = [];
model.zrange = [];
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
cdata = permute(cdata, [2 1 3 4]);
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


opts = {'Parent',ax,'CDataMapping',[],'AlphaDataMapping',[],'FaceColor','texturemap','FaceAlpha','texturemap','EdgeColor','none','tag',model.tag};

if ndims(cdata) > 3
    opts{4} = 'direct';
else
    if isempty(model.cscale)
       %need to scale manually to use cdata as alpha data too !
       cdata = cdata - min(cdata(:));
       cdata = cdata / max(cdata(:));
       caxis(ax, [0, 1]); 
       %caxis(ax, [min(cdata(:)), max(cdata(:))])
    else
       cdata = cdata - model.cscale(1);
       cdata(cdata < 0) = 0;
       cdata = cdata / (model.cscale(2) - model.cscale(1));
       cdata(cdata > 1) = 1;
       caxis(ax, model.cscale)
    end
    opts{4} = 'scaled';
end


if isempty(model.alpha)
    alpha = cdata;
    if ndims(alpha) > 3
        alpha = sqrt(sum(double(alpha).^2, 4));
        alpha = alpha - min(alpha(:));
        alpha = alpha / max(alpha(:));
    end
    opts{6} = 'none';
    
    %max(alpha(:))
    %min(alpha(:))
else
   alpha = model.alpha;
   
   if ndims(alpha) == 3
      alpha = permute(alpha, [2 1 3]);
   else
      alpha = permute(cdata, [2 1 3 4]);
      alpha = sqrt(sum(double(alpha).^2, 4));
      alpha = alpha - min(alpha(:));
      alpha = alpha / max(alpha(:));
   end
   
   if ~isequal(csiz(1:3), size(alpha))
      error('implot3d: Incorrect size of alpha!');
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
   x = [model.xrange(1), model.xrange(2); model.xrange(1), model.xrange(2)];
   y = [model.yrange(1), model.yrange(1); model.yrange(2), model.yrange(2)];
   z = [model.zrange(1), model.zrange(1); model.zrange(1), model.zrange(1)];
   diff = model.zrange(2)-model.zrange(1);
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
   x = [model.xrange(1), model.xrange(1); model.xrange(1), model.xrange(1)];
   y = [model.yrange(1), model.yrange(1); model.yrange(2), model.yrange(2)];
   z = [model.zrange(1), model.zrange(2); model.zrange(1), model.zrange(2)];
   diff = model.xrange(2)-model.xrange(1);
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
   x = [model.xrange(1), model.xrange(1); model.xrange(2), model.xrange(2)];
   y = [model.yrange(1), model.yrange(1); model.yrange(1), model.yrange(1)];
   z = [model.zrange(1), model.zrange(2); model.zrange(1), model.zrange(2)];
   diff = model.yrange(2)-model.yrange(1);
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
