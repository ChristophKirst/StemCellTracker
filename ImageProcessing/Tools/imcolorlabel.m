function h = imcolorlabel(varargin)
%
% imcolorlabel
% imcolorlabel off
% imcolorlabel(label, colors)
% imcolorlabel(haxes, ...)
% h = imcolorlabel(...)
% 
% description:
%     draws color label onto the images (e.g use to label channel colors)
% 
% inputs:
%     paramter struct with entries:
%     haxes             handle to the axes (defaults to current axes)
%     location          location of the scalebar: 'ne', 'nw', 'se', 'sw' or [x,y] ('ne' default)
%     color             color of the different label as rgb values in a cell
%     fontweight        font weight ('normal')
%     textalignment     alignment of the text e.g. {'center', 'top'} for horizontal bars ([] = automatic)
%     separator         separator for the labels ('/')
%
% output:
%     h                 reference to scalebar object
%
% See also: imscalebar

narg = nargin;
if narg > 0
   if ishandle(varargin{1}) %haxes
      varargin{length(varargin)+1}='haxes';
      varargin{length(varargin)+1}=varargin{1};
      varargin = varargin(2:end);
      narg = narg -1;
   end
   
   % delete label if requested
   if narg == 1 && strcmpi(varargin{1},'off')
      cb = findobj(gca,'tag','colorlabel');
      if ~isempty(cb)
         delete(cb);
      end
      return;
   end
end


if narg > 0 && iscell(varargin{1})
   varargin{length(varargin)+1}='label';
   varargin{length(varargin)+1}=varargin{1};
   varargin = varargin(2:end);
   narg = narg -1;
end

if narg > 0 && iscell(varargin{1})
   varargin{length(varargin)+1}='color';
   varargin{length(varargin)+1}=varargin{1};
   varargin = varargin(2:end);
   %narg = narg - 1;
end


% params
param = parseParameter(varargin);

hAxes            = getParameter(param, 'haxes', gca);
location         = getParameter(param, 'location', 'ne');
textalgn         = getParameter(param, 'textalign', []); 
color           = getParameter(param, 'color', {[1,1,1]});
fontweight       = getParameter(param, 'fontweight', 'normal');
label            = getParameter(param, 'label', '');
separator        = getParameter(param, 'separator', ' ');

if isempty(label)
   return
end
if ~iscell(label)
   label = {label};
end
if isempty(color)
   color = {[1,1,1]};
end
if ~iscell(color)
   color = {color};
end


cb = findobj(hAxes,'tag','colorlabel');
if ~isempty(cb)
   delete(cb);
end


if ~ishandle(hAxes)
   error 'haxes is not a valid Axes handle' ;
elseif ~strcmpi(get(hAxes,'type'),'axes')
   error 'haxes is not a valid Axes handle';
elseif strcmpi(get(hAxes,'tag'),'colorbar')
   error 'haxes is a handle to a colorbar';
end
          
directions = {'ne', 'nw', 'se', 'sw'};

if ~(numel(location)==2 && isnumeric(location)) && ...
      isempty(strcmpi(location,directions))
   error 'unrecognised value for location'
end

% CHECK IF VIEW IS IN 2D
[~, el] = view(hAxes);
if el~=90
   error 'The Axes must be in 2D view'
end
    
%GET IMAGE AND AXES DATA
axeslims = [get(hAxes,'xlim')' get(hAxes,'ylim')'];
axesdir = [1 1];
if strcmpi(get(hAxes,'XDir'),'reverse')
   axeslims(:,1) = flipud(axeslims(:,1));
   axesdir(1) = -1;
end
if strcmpi(get(hAxes,'YDir'),'reverse')
   axeslims(:,2) = flipud(axeslims(:,2));
   axesdir(2) = -1;
end

rng1 = range2(axeslims(:,1));
rng2 = range2(axeslims(:,2));


%SET UP POSITIONING
fac = 0.05;

if ischar(location)
   switch location
      case 'ne'
         anchor = [axeslims(2,1) - axesdir(1) * rng1 * fac, ...
                   axeslims(2,2) - axesdir(2) * rng2 * fac];
      case 'nw'
         anchor = [axeslims(1,1) + axesdir(1) * rng1 * fac, ...
                   axeslims(2,2) - axesdir(2) * rng2 * fac];
      case 'sw'
         anchor = [axeslims(1,1) + axesdir(1) * rng1 * fac, ...
                   axeslims(1,2) + axesdir(2) * rng2 * fac];
      case 'se'         
         anchor = [axeslims(2,1) - axesdir(1) * rng1 * fac, ...
                   axeslims(1,2) + axesdir(2) * rng2 * fac];
   end
else
   anchor = location;
end

if isempty(textalgn)
   switch location
      case 'ne'
         textalgn = {'top', 'right'};
      case 'nw'
         textalgn = {'top', 'left'};
      case 'se'
         textalgn = {'bottom', 'right'};
      case 'sw'         
         textalgn = {'bottom', 'left'};
      otherwise
         textalgn = {'bottom', 'left'};
   end
end

% GENERATE COLORED TEXT

txt = '';
for i = 1:length(label)
   txt = [txt, sprintf('{\\color[rgb]%s %s}', var2char(num2cell(color{i})), label{i})]; %#ok<AGROW>
   if i < length(label)
      txt = [txt, separator]; %#ok<AGROW>
   end
end
      


% DRAW SCALEBAR
set(hAxes,'xlimmode','manual','ylimmode','manual');
hg = hggroup('tag','colorlabel');

text(anchor(1), anchor(2), 0, txt, 'verticalalignment',textalgn{1},'horizontalalignment',textalgn{2}, 'fontweight', fontweight, 'parent', hg);

if nargout>0
   h = hg;
end
    
% SETUP DELETE CALLBACK
set(hg,'DeleteFcn',@deleteColorLabel)

% SETUP LISTENER TO RESET SCALEBAR ON CHANGE OF AXES LIMITS
hL(1) = addlistener(hAxes,'YLim','PostSet',@(src,event) resetColorLabel(src,event,hg));

% SET USERDATA

if isfield(param, 'haxes')
   param = rmfield(param, 'haxes');
end

udata = get(hg,'UserData');
udata = parseParameter(udata, 'CLAnchorRatio',[(anchor(1)-min(axeslims(:,1)))/range2(axeslims(:,1)) (anchor(2)-min(axeslims(:,2)))/range2(axeslims(:,2))],...
         'CLListeners',hL, 'CLParameter', param);
set(hg,'UserData',udata);
     
end

% CALLBACK FUNCTIONS
function deleteColorLabel(src,event) %#ok<INUSD>
   udata = get(src,'UserData');
   if isentry(udata, 'CBListeners')
      l = udata.('CBListeners');
      if ~isempty(l)
         delete(l);
      end
   end
end

function resetColorLabel(src,event,SB) %#ok<INUSL>
   udata = get(SB,'UserData');
   hAxes = get(SB,'parent');

   delete(SB);

   axeslims = [get(hAxes,'xlim')' get(hAxes,'ylim')'];

   anchorratio = udata.('CBAnchorRatio');
   location = [anchorratio(1)*range2(axeslims(:,1))+axeslims(1,1) anchorratio(2)*range2(axeslims(:,2))+axeslims(1,2)];
   
   param = udata.('Parameter');

   imcolorlabel(hAxes, param, 'Location',location);
end

function x = range2(x)
   x = max(x) - min(x);
end
          
        