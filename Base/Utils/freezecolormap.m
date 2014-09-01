function freezecolormap(varargin)
%
% freezecolormap()
% freezecolormap(axh)
%
% description:
%   freeze colormap in a subplot
%
% input:
%   axh   (optional) the axis handle for which to freeze the color map
%
% See also: implot
%

% Based on freezeColors by, John Iversen, 2005
% Modifications: storing color data for large/many images is not usefull, so we skip the option to undo freezing

[h, nancolor] = checkArgs(varargin);

%gather all children with scaled or indexed CData
cdatah = getCDataHandles(h);

%current colormap
cmap = colormap;
nColors = size(cmap,1);
cax = caxis;

% convert object color indexes into colormap to true-color data using 
%  current colormap
for hh = cdatah',
    g = get(hh);
    
    %preserve parent axis clim
    parentAx = getParentAxes(hh);
    originalClim = get(parentAx, 'clim');    
   
    %   Note: Special handling of patches: For some reason, setting
    %   cdata on patches created by bar() yields an error,
    %   so instead we'll set facevertexcdata instead for patches.
    if ~strcmp(g.Type,'patch'),
        cdata = g.CData;
    else
        cdata = g.FaceVertexCData; 
    end
    
    %get cdata mapping (most objects (except scattergroup) have it)
    if isfield(g,'CDataMapping'),
        scalemode = g.CDataMapping;
    else
        scalemode = 'scaled';
    end
    
    %save original indexed data for use with unfreezeColors
    siz = size(cdata);
    %we dont weant to store all the data !
    %setappdata(hh, appdatacode, {cdata scalemode});

    %convert cdata to indexes into colormap
    if strcmp(scalemode,'scaled'),
        %4/19/06 JRI, Accommodate scaled display of integer cdata:
        %       in MATLAB, uint * double = uint, so must coerce cdata to double
        %       Thanks to O Yamashita for pointing this need out
        idx = ceil( (double(cdata) - cax(1)) / (cax(2)-cax(1)) * nColors);
    else %direct mapping
        idx = cdata;
        %10/8/09 in case direct data is non-int (e.g. image;freezeColors)
        % (Floor mimics how matlab converts data into colormap index.)
        % Thanks to D Armyr for the catch
        idx = floor(idx);
    end
    
    %clamp to [1, nColors]
    idx(idx<1) = 1;
    idx(idx>nColors) = nColors;

    %handle nans in idx
    nanmask = isnan(idx);
    idx(nanmask)=1; %temporarily replace w/ a valid colormap index

    %make true-color data--using current colormap
    realcolor = zeros(siz);
    for i = 1:3,
        c = cmap(idx,i);
        c = reshape(c,siz);
        c(nanmask) = nancolor(i); %restore Nan (or nancolor if specified)
        realcolor(:,:,i) = c;
    end
    
    %apply new true-color color data
    
    %true-color is not supported in painters renderer, so switch out of that
    if strcmp(get(gcf,'renderer'), 'painters'),
        set(gcf,'renderer','zbuffer');
    end
    
    %replace original CData with true-color data
    if ~strcmp(g.Type,'patch'),
        set(hh,'CData',realcolor);
    else
        set(hh,'faceVertexCData',permute(realcolor,[1 3 2]))
    end
    
    %restore clim (so colorbar will show correct limits)
    if ~isempty(parentAx),
        set(parentAx,'clim',originalClim)
    end
    
end %loop on indexed-color objects


% ============================================================================ %
% Local functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getCDataHandles -- get handles of all descendents with indexed CData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hout = getCDataHandles(h)
% getCDataHandles  Find all objects with indexed CData

%recursively descend object tree, finding objects with indexed CData
% An exception: don't include children of objects that themselves have CData:
%   for example, scattergroups are non-standard hggroups, with CData. Changing
%   such a group's CData automatically changes the CData of its children, 
%   (as well as the children's handles), so there's no need to act on them.

narginchk(1,1)

hout = [];
if isempty(h),return;end

ch = get(h,'children');
for hh = ch'
    g = get(hh);
    if isfield(g,'CData'),     %does object have CData?
        %is it indexed/scaled?
        if ~isempty(g.CData) && isnumeric(g.CData) && size(g.CData,3)==1, 
            hout = [hout; hh]; %#ok<AGROW> %yes, add to list
        end
    else %no CData, see if object has any interesting children
            hout = [hout; getCDataHandles(hh)]; %#ok<AGROW>
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getParentAxes -- return handle of axes object to which a given object belongs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hAx = getParentAxes(h)
% getParentAxes  Return enclosing axes of a given object (could be self)

narginchk(1,1)
%object itself may be an axis
if strcmp(get(h,'type'),'axes'),
    hAx = h;
    return
end

parent = get(h,'parent');
if (strcmp(get(parent,'type'), 'axes')),
    hAx = parent;
else
    hAx = getParentAxes(parent);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% checkArgs -- Validate input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h, nancolor] = checkArgs(args)
% checkArgs  Validate input arguments to freezeColors

nargs = length(args);
if nargs < 0 || nargs > 3
   error('wrong number of argumetns!')
end

%grab handle from first argument if we have an odd number of arguments
if mod(nargs,2),
    h = args{1};
    if ~ishandle(h),
        error('JRI:freezeColors:checkArgs:invalidHandle',...
            'The first argument must be a valid graphics handle (to an axis)')
    end
    % 4/2010 check if object to be frozen is a colorbar
    if strcmp(get(h,'Tag'),'Colorbar'),
      if ~exist('cbfreeze.m', 'file'),
        warning('JRI:freezeColors:checkArgs:cannotFreezeColorbar',...
            ['You seem to be attempting to freeze a colorbar. This no longer'...
            'works. Please read the help for freezeColors for the solution.'])
      else
        cbfreeze(h);
        return
      end
    end
    args{1} = [];
    nargs = nargs-1;
else
    h = gca;
end

%set nancolor if that option was specified
nancolor = [nan nan nan];
if nargs == 2,
    if strcmpi(args{end-1},'nancolor'),
        nancolor = args{end};
        if ~all(size(nancolor)==[1 3]),
            error('JRI:freezeColors:checkArgs:badColorArgument',...
                'nancolor must be [r g b] vector');
        end
        nancolor(nancolor>1) = 1; nancolor(nancolor<0) = 0;
    else
        error('JRI:freezeColors:checkArgs:unrecognizedOption',...
            'Unrecognized option (%s). Only ''nancolor'' is valid.',args{end-1})
    end
end