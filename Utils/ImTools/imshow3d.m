function [model] = imshow3d(varargin)
%
% [model] = imshow3d(varargin)
%
% description:
%    im3d uses the orthogonal plane 2-D texture mapping technique for 
%    volume rending 3-D data in OpenGL. Use the 'texture' option to fine 
%    tune the texture mapping technique. This function is best used with
%    fast OpenGL hardware.
%
% usage:
% imshow3d                 Provide a demo of functionality.
%
% H = imshow3d('CData',data)   Create volume render object from input 
%                              3-D data. Use interp3 on data to increase volume
%                              rendering resolution. Returns a struct 
%                              encapsulating the pseudo-volume rendering object.
%                              XxYxZ array represents scaled colormap indices.
%                              XxYxZx3 array represents truecolor RGB values for
%                              each voxel (along the 4th dimension).
%
% imshow3d(...,'Alpha',alpha)  XxYxZ array of alpha values for each voxel, in
%                              range [0,1]. Default: data (interpreted as
%                              scaled alphamap indices).
%
% imshow3d(...,'Parent',axH)   Specify parent axes. Default: gca.
%
% imshow3d(...,'HRange',y)     1x2 h-axis bounds. Default: [1 size(data, 1)].
% imshow3d(...,'WRange',x)     1x2 w-axis bounds. Default: [1 size(data, 2)].
% imshow3d(...,'LRange',z)     1x2 z-axis bounds. Default: [1 size(data, 3)].
% imshow3d(...,'Range', r)     1x6 h,w,z bounds [xmin xmax; xmin, ymax;  zmin, zmax];
% imshow3d(...,'BoxRatios', r) 1x3 ratio array: coordinate ranges are [1 r(1)*size(data,1)], [1 r(2)*size(data,1)]...
%
% imshow3d(...,'texture','2D') Only render texture planes parallel to nearest
%                              orthogonal viewing plane. Requires doing
%                              imshow3d(h) to refresh if the view is rotated
%                              (i.e. using cameratoolbar).
%
% imshow3d(...,'texture','3D') Default. Render x,y,z texture planes
%                              simultaneously. This avoids the need to
%                              refresh the view but requires faster OpenGL
%                              hardware peformance.
%
% imshow3d(H)                  Refresh view. Updates rendering of texture planes 
%                              to reduce visual aliasing when using the 'texture'='2D' option.
%
% note:
% Use interp3 on input date to increase/decrease resolution of data
%
% See also alphamap, colormap, opengl, isosurface

% Copyright Joe Conti, 2004
% modified by C. Kirst 2014 to use pixel coordinates, 
%                           correct for similar appearance as in montage, and using data ranges and box ratios 

if isstruct(varargin{1})
    model = varargin{1};
    if length(varargin) > 1
       varargin = {varargin{2:end}}; %#ok<CCAT1>
    end
else
    model = localGetDefaultModel;
end

if ~ischar(varargin{1})
   varargin = {'cdata', varargin{:}};
end

if length(varargin)>1
  for n = 1:2:length(varargin)
    switch(lower(varargin{n}))
        case 'cdata'
           model.cdata = varargin{n+1};
        case 'parent'
           model.parent = varargin{n+1};
        case 'texture'
           model.texture = varargin{n+1};
        case 'alpha'
           model.alpha = varargin{n+1};
        case 'xrange'
           model.xdata = varargin{n+1}([1 end]);
        case 'yrange'
           model.ydata = varargin{n+1}([1 end]);
       case 'zrange'
           model.zdata = varargin{n+1}([1 end]);
       case 'hrange'
           model.xdata = varargin{n+1}([1 end]);
       case 'wrange'
           model.ydata = varargin{n+1}([1 end]);
       case 'lrange'
           model.zdata = varargin{n+1}([1 end]);
        case 'range'  % assume pixel coordinates
           r = varargin{n+1}(:);
           model.xdata = r(1:2);         
           model.ydata = r(3:4);         
           model.zdata = r(5:6);
    end
  end
end

% Define [x,y,z]data
siz = size(model.cdata);
if isempty(model.xdata)
    model.xdata = [0 siz(1)] + 0.5;  % we use pixel coordinates,model here stores spatial coordinates [h,w] = [y,x]
end
if isempty(model.ydata)
    model.ydata = [0 siz(2)] + 0.5;
end
if isempty(model.zdata)
    model.zdata = [0 siz(3)] + 0.5;
end


if length(varargin)>1
  for n = 1:2:length(varargin)
    switch(lower(varargin{n}))
        case 'boxratios'
           r = varargin{n+1}(:);
           model.xdata = [model.xdata(1) r(1) * (model.xdata(2)-model.xdata(1)) + model.xdata(1)];
           model.ydata = [model.ydata(1) r(2) * (model.ydata(2)-model.ydata(1)) + model.ydata(1)];   
           model.zdata = [model.zdata(1) r(3) * (model.zdata(2)-model.zdata(1)) + model.zdata(1)];
    end    
  end
end

if isempty(model.parent)
    model.parent = gca;
end

[model] = local_draw(model);

set(gca,'YDir','Reverse');
xlabel('w'); ylabel('h'); zlabel('l')

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
model.texture = '3D';
tag = tempname;
model.tag = ['vol3d_' tag(end-11:end)];

end

%------------------------------------------%
function [model,ax] = local_draw(model)

cdata = model.cdata;
cdata = permute(cdata, [2 1 3]);
siz = size(cdata);

try
   delete(model.handles);
catch
end

ax = model.parent;
cam_dir = camtarget(ax) - campos(ax);
[~,ind] = max(abs(cam_dir));

opts = {'Parent',ax,'cdatamapping',[],'alphadatamapping',[],'facecolor','texturemap','edgealpha',0,'facealpha','texturemap','tag',model.tag};

if ndims(cdata) > 3
    opts{4} = 'direct';
else
    cdata = double(cdata);
    opts{4} = 'scaled';
end

if isempty(model.alpha)
    alpha = cdata;
    if ndims(model.cdata) > 3
        alpha = sqrt(sum(double(alpha).^2, 4));
        alpha = alpha - min(alpha(:));
        alpha = 1 - alpha / max(alpha(:));
    end
    opts{6} = 'scaled';
else
    alpha = model.alpha;
    if ~isequal(siz(1:3), size(alpha))
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

is3DTexture = strcmpi(model.texture,'3D');
handle_ind = 1;

% Create z-slice
if(ind==3 || is3DTexture )    
  x = [model.xdata(1), model.xdata(2); model.xdata(1), model.xdata(2)];
  y = [model.ydata(1), model.ydata(1); model.ydata(2), model.ydata(2)];
  z = [model.zdata(1), model.zdata(1); model.zdata(1), model.zdata(1)];
  diff = model.zdata(2)-model.zdata(1);
  delta = diff/size(cdata,3);
  for n = 1:size(cdata,3)

   cslice = squeeze(cdata(:,:,n,:));
   aslice = double(squeeze(alpha(:,:,n)));
   h(handle_ind) = surface(x,y,z,cslice,'alphadata',aslice,opts{:});
   z = z + delta;
   handle_ind = handle_ind + 1;
  end

end

% Create x-slice
if (ind==1 || is3DTexture ) 
  x = [model.xdata(1), model.xdata(1); model.xdata(1), model.xdata(1)];
  y = [model.ydata(1), model.ydata(1); model.ydata(2), model.ydata(2)];
  z = [model.zdata(1), model.zdata(2); model.zdata(1), model.zdata(2)];
  diff = model.xdata(2)-model.xdata(1);
  delta = diff/size(cdata,2);
  for n = 1:size(cdata,2)

   cslice = squeeze(cdata(:,n,:,:));
   aslice = double(squeeze(alpha(:,n,:)));
   h(handle_ind) = surface(x,y,z,cslice,'alphadata',aslice,opts{:});
   x = x + delta;
   handle_ind = handle_ind + 1;
  end
end

% Create y-slice
if (ind==2 || is3DTexture)
  x = [model.xdata(1), model.xdata(1); model.xdata(2), model.xdata(2)];
  y = [model.ydata(1), model.ydata(1); model.ydata(1), model.ydata(1)];
  z = [model.zdata(1), model.zdata(2); model.zdata(1), model.zdata(2)];
  diff = model.ydata(2)-model.ydata(1);
  delta = diff/size(cdata,1);
  for n = 1:size(cdata,1)

   cslice = squeeze(cdata(n,:,:,:));
   aslice = double(squeeze(alpha(n,:,:)));
   h(handle_ind) = surface(x,y,z,cslice,'alphadata',aslice,opts{:});
   y = y + delta;
   handle_ind = handle_ind + 1;
  end
end

model.handles = h;

end


% function demo_imshow3d
% figure;
% load mri.mat
% vol3d('cdata', squeeze(D), 'xdata', [0 1], 'ydata', [0 1], 'zdata', [0 0.7]);
% colormap(bone(256));
% alphamap([0 linspace(0.1, 0, 255)]);
% axis equal off
% set(gcf, 'color', 'w');
% view(3);
% 
% end