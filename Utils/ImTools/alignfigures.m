function alignfigures(figs, style, screenarea)
%
% tilefigures(figs)
%
% description:
%    tiles / aranges figures on screen for better clarity
%
% input:
%    figs        specific figures (all)
%    style       how to align 
%    screenarea  the screen area to arage figures in
%

% default values
default_style = 'tile';

screen_top = 837;  % set this in multiple monitor setup if dection is bad
screen_height = 736;  % set this in multiple monitor setup if dection is bad
%screen_width = 837;  % set this in multiple monitor setup if dection is bad
screen_left = 0;

if nargin < 1
   figs = 'all';
end
if nargin < 2
   if ischar(figs)
      switch figs
         case 'tile'
            figs = 'all';
            style = 'tile';
         case 'offset'
            figs = 'all';
            style = 'offset';
         otherwise
            figs = 'all';
            style = default_style;
      end
   else
      style = default_style;
   end
end
 
if isequal(figs, 'all')
   figs = findobj('Type', 'figure');
   figs	= sort(figs);
end

if isempty(figs)
   disp('tilefigures: no figures to arrange!');
   return
end



if nargin < 3
   screenarea = 'monitor 1';
end

units = 'pixels';
old_root_units = get(0, 'Units');	
set(0, 'Units', units);	

if ischar(screenarea)
   switch screenarea
      case 'all'
         screenarea = get(0, 'ScreenSize');
      case 'monitor 1'
         mp = get(0,'MonitorPositions');
         screenarea = mp(1,:);  
      case 'monitor 2'
         mp = get(0,'MonitorPositions');
         screenarea = mp(2,:); 
      otherwise 
         error('alignfigures: screen area specification not supported!')
   end
end

set(0, 'Units', old_root_units);
sx = screenarea(1);
sy = screenarea(2);
sw = screenarea(3);			
sh = screenarea(4);

% multi monitors
if exist('screen_left', 'var')
   sx = screen_left;
end
if exist('screen_top', 'var')
   sy = screen_top;
end
if exist('screen_height', 'var')
   sh = screen_height;
end
if exist('screen_width', 'var')
   sw = screen_width;
end



nfigs = length(figs);

switch style
   case 'offset'
      
      fh = figs(end);
      old_fig_units = get(fh, 'Units');
	   set(fh, 'Units', units);
	   fp = get(fh, 'OuterPosition');
      set(fh, 'Units', old_fig_units);
      
      off = (sw - fp(3))/nfigs;
      off = min(off, (sh - fp(4))/nfigs);
      
      top = sy + sh;
      
      for i = 1 : nfigs
         fh = figs(i);
         old_fig_units = get(fh, 'Units');
         set(fh, 'Units', units);
         fp = get(fh, 'OuterPosition');
         fp = [sx + (i-1)*off top - (i-1)*off fp(3) fp(4)];
         set(fh, 'OuterPosition', fp);
	      set(fh, 'Units', old_fig_units);
	      figure(fh); % raise figure.
      end  
      
   case 'tile'
          
      nh = ceil(sqrt(nfigs));
      nv = ceil(nfigs / nh);

      % some spacings
      hspc   = 0;  % horisontal space.
      topspc = 0;  % space above top figure.
      medspc = 10; % space between figures.
      botspc = 38; % pace below bottom figure.

      fw = (sw - (nh + 1) * hspc) / nh;
      fh = (sh - (topspc + botspc) - (nv - 1) * medspc) / nv;

      for row = 1 : nv
         for col = 1 : nh
            idx = (row - 1) * nh + col;
            if idx <= nfigs
                fx =  (col - 1) * fw + col * hspc;
                fy = sy + sh - topspc - (row-1) * fh + (row - 1) * medspc;

                fp = [ fx fy fw fh ];
                fi = figs(idx);
                old_fig_units = get(fi, 'Units');
                set(fi, 'Units', units);
                set(fi, 'OuterPosition', fp);
                set(fi, 'Units', old_fig_units);
                figure(fi);
            end
         end
      end
 
   otherwise
      error('alignfigures: style specification not supported!')
end

end

