function h = imscalebar(varargin)
%
% imscalebar
% imscalebar off
% imscalebar(param)
% imscalebar(haxes, param)
% h = imscalebar(...)
% 
% description:
%    draws a scalebar on the axes and returns handle to the scalebar using the specified parameter
% 
% inputs:
%     paramter struct with entries:
%     haxes             handle to the axes (defaults to current axes)
%     length            length of scale bar in final units (1)
%     scale             scale of the image, i.e. 1 unit in image  = scale units in real (1)
%     location          location of the scalebar: 'ne', 'nw', 'se', 'sw' or [x,y] ('ne' default)
%     orientation       orientation 'horizontal' or 'vertical' ('horizontal')
%     color             color of scalebar (default is [0 0 0])
%     fontweight        font weight ('normal')
%     linewidth         axesdir(2)*0.1*scalelengthwidth of scale bar line (0.5)
%     unit              string containing unit e.g. 'mm'
%     textoffset        offset of text from center of scale bar in percent of image  ([] = automatic)
%     textalignment     alignment of the text e.g. {'center', 'top'} for horizontal bars ([] = automatic)
%     edgewidth         edage width in percent of the image (0.01)
%
% output:
%     h                 reference to scalebar object
%
% Note: imscalebar sets the XLimMode and YLimMode of the axes to manual.

% Acknowledgment: Based on scalebar by Amanda Ng

narg = nargin;
if narg > 0
   if ishandle(varargin{1}) %haxes
      varargin{length(varargin)+1}='haxes';
      varargin{length(varargin)+1}=varargin{1};
      varargin = varargin(2:end);
      narg = narg -1;
   end
   
   % delete scale bar if requested
   if narg == 1 && strcmpi(varargin{1},'off')
      sb = findobj(gca,'tag','scalebar');
      if ~isempty(sb)
         delete(sb);
      end
      return;
   end
end

% params
param = parseParameter(varargin);

hAxes            = getParameter(param, 'haxes', gca);
scale            = getParameter(param, 'scale', 1);
len              = getParameter(param, 'length', 1);
location         = getParameter(param, 'location', 'ne');
orientation      = getParameter(param, 'orientation', 'horizontal');
textoffset       = getParameter(param, 'textoffset', []); 
textalgn         = getParameter(param, 'textalign', []); 
colour           = getParameter(param, 'color', [1,1,1]);
linewidth        = getParameter(param, 'linewidth', 0.5);
fontweight       = getParameter(param, 'fontweight', 'normal');
unitstring       = getParameter(param, 'unit', '');
edgewidth        = getParameter(param, 'edgewidth', 0.01);


sb = findobj(hAxes,'tag','scalebar');
if ~isempty(sb)
   delete(sb);
end


% constants
directions = {'ne', 'nw', 'se', 'sw'};

if ~ishandle(hAxes)
   error 'haxes is not a valid Axes handle' ;
elseif ~strcmpi(get(hAxes,'type'),'axes')
   error 'haxes is not a valid Axes handle';
elseif strcmpi(get(hAxes,'tag'),'colorbar')
   error 'haxes is a handle to a colorbar';
end
                   
if ~isnumeric(scale)
   error 'scale must be a numeric value';
end

if ~isnumeric(len)
   error 'length must be a numeric value'
end

if ~(numel(location)==2 && isnumeric(location)) && ...
      isempty(strcmpi(location,directions))
   error 'unrecognised value for location'
end

if numel(colour)~=3 || ~isnumeric(colour)
   error 'colour must be a 1x3 representation of an RGB colour'
end

if ~ischar(unitstring) || isempty(unitstring)
   error('''Unit'' must be followed by a string')
end
unitstring = [' ' strtrim(unitstring)];

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

%COMPUTE SCALE
scalelength = len/scale;

%SET UP POSITIONING
fac = 0.05;

if ischar(location)
   switch location
      case 'ne'
         anchor = [axeslims(2,1) - axesdir(1) * rng1 * fac, ...
                   axeslims(2,2) - axesdir(2) * rng2 * fac];
         if strcmp(orientation, 'horizontal') 
            sbdir  = [-1,0];
         else
            sbdir  = [0,-1];
         end
      case 'nw'
         anchor = [axeslims(1,1) + axesdir(1) * rng1 * fac, ...
                   axeslims(2,2) - axesdir(2) * rng2 * fac];
         if strcmp(orientation, 'horizontal') 
            sbdir  = [1,0];
         else
            sbdir  = [0,-1];
         end
      case 'sw'
         anchor = [axeslims(1,1) + axesdir(1) * rng1 * fac, ...
                   axeslims(1,2) + axesdir(2) * rng2 * fac];
         if strcmp(orientation, 'horizontal') 
            sbdir  = [1,0];
         else
            sbdir  = [0,1];
         end
                
      case 'se'         
         anchor = [axeslims(2,1) - axesdir(1) * rng1 * fac, ...
                   axeslims(1,2) + axesdir(2) * rng2 * fac];
         if strcmp(orientation, 'horizontal') 
            sbdir  = [-1,0];
         else
            sbdir  = [0,1];
         end
   end
else
   anchor = location;
   if strcmp(orientation, 'horizontal')
      sbdir  = [1,0];
   else
      sbdir  = [0,1];
   end
end

sbdir = axesdir .* sbdir;

linepos = [anchor; anchor + sbdir * scalelength];

if strcmp(orientation, 'horizontal')
   edir  = [0,1] .* edgewidth .* [rng1, rng2];
else
   edir  = [1,0] .* edgewidth .* [rng1, rng2];
end
edir = [edir; -edir];

if isempty(textoffset)
   if strcmp(orientation, 'horizontal')
      textoffset  = [0,-1] * rng2 * 0.01;
   else
      textoffset  = [1, 0] * rng2 * 0.01;
   end
end

if isscalar(textoffset)
   if strcmp(orientation, 'horizontal')
      textoffset  = [0,textoffset] * rng2;
   else
      textoffset  = [textoffset, 0] * rng2;
   end
end

if isempty(textalgn)
   if strcmp(orientation, 'horizontal')
      if textoffset(2) > 0
         textalgn = {'bottom', 'center'};
      else
         textalgn = {'top', 'center'};
      end
   else
      if textoffset(1) > 0
         textalgn = {'center', 'left'};
      else
         textalgn = {'center', 'right'};
      end
   end
end


            
% DRAW SCALEBAR
set(hAxes,'xlimmode','manual','ylimmode','manual');
hg = hggroup('tag','scalebar');
line(linepos(:,1), linepos(:,2), 'color', colour, 'linewidth', linewidth, 'parent', hg);
line(linepos(1,1) + edir(:,1), linepos(1,2) + edir(:,2), 'color', colour, 'linewidth', linewidth, 'parent', hg);
line(linepos(2,1) + edir(:,1), linepos(2,2) + edir(:,2), 'color', colour, 'linewidth', linewidth, 'parent', hg);
%text(linepos(1,1),linepos(1,2),0,'0','verticalalignment',textalignment{1},'horizontalalignment',textalignment{2}, 'color', colour, 'fontweight', fontweight, 'parent', hg);
text(linepos(1,1) + (linepos(2,1)-linepos(1,1))/2 + textoffset(1), linepos(1,2) + (linepos(2,2)-linepos(1,2))/2 + textoffset(2), 0, [num2str(len) unitstring], ...
     'verticalalignment',textalgn{1},'horizontalalignment',textalgn{2}, ...
     'color', colour, 'fontweight', fontweight, 'parent', hg);

if nargout>0
   h = hg;
end
    
% SETUP DELETE CALLBACK
set(hg,'DeleteFcn',@deleteScaleBar)

% SETUP LISTENER TO RESET SCALEBAR ON CHANGE OF AXES LIMITS
hL(1) = addlistener(hAxes,'YLim','PostSet',@(src,event) resetScaleBar(src,event,hg));

% SET USERDATA

if isfield(param, 'haxes')
   param = rmfield(param, 'haxes');
end

udata = get(hg,'UserData');
udata = parseParameter(udata, 'SBAnchorRatio',[(anchor(1)-min(axeslims(:,1)))/range2(axeslims(:,1)) (anchor(2)-min(axeslims(:,2)))/range2(axeslims(:,2))],...
         'SBListeners',hL, 'SBParameter', param);
set(hg,'UserData',udata);
     
end

% CALLBACK FUNCTIONS
function deleteScaleBar(src,event) %#ok<INUSD>
   udata = get(src,'UserData');
   if isentry(udata, 'SBListeners')
      l = udata.('SBListeners');
      if ~isempty(l)
         delete(l);
      end
   end
end

function resetScaleBar(src,event,SB) %#ok<INUSL>
   udata = get(SB,'UserData');
   hAxes = get(SB,'parent');

   delete(SB);

   axeslims = [get(hAxes,'xlim')' get(hAxes,'ylim')'];

   anchorratio = udata.('SBAnchorRatio');
   location = [anchorratio(1)*range2(axeslims(:,1))+axeslims(1,1) anchorratio(2)*range2(axeslims(:,2))+axeslims(1,2)];
   
   param = udata.('SBParameter');

   imscalebar(hAxes, param, 'Location',location);
end

function x = range2(x)
   x = max(x) - min(x);
end
          
        