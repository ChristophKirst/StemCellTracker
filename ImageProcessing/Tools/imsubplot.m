function h_subplot = imsubplot(n, m, p, varargin)
%
% h_subplot = imsubplot(n, m, p, param)
%
% description:
%    creates subplot with small spacings
%
% input:
%     n,m,p     as in subplot(m,n,p), n,m are revers3ed for pql consistency and numbering starts left low to right up
%     param     parameter struct with entries
%               gap       (optional) spacing between plots
%               marg      (optional) border spacings
%               subplot   other inputs to subplot {}
%
% output:
%     h_subplot handle to subplot
%
% See also: subplot

param = parseParameter(varargin{:});

gap = getParameter(param, 'gap', 0.005);
marg= getParameter(param,'marg', 0.005);
varg= getParameter(param,'subplot', {});

if numel(marg)==2
    marg_w = [marg(1) marg(1)];
    marg_h = [marg(2) marg(2)];
elseif numel(marg) == 4
    marg_w = [marg(1) marg(2)];
    marg_h = [marg(3) marg(4)];
else
    marg_w = [marg(1) marg(1)];
    marg_h = [marg(1) marg(1)];
end
    
if numel(gap)==1; 
    gap_w = gap;
    gap_h = gap;
else
    gap_w = gap(1);
    gap_h = gap(2);
end

[subplot_col, subplot_row] = ind2sub([n,m],p);
subplot_row = m - subplot_row + 1;

height = (1 - marg_h(1) - marg_h(2) - (m-1) * gap_h ) / m; % single subplot height
width  = (1 - marg_w(1) - marg_w(2) - (n-1) * gap_w ) / n; % single subplot width

% note subplot suppors vector p inputs- so a merged subplot of higher dimensions will be created
subplot_cols = 1+max(subplot_col)-min(subplot_col); % number of column elements in merged subplot 
subplot_rows = 1+max(subplot_row)-min(subplot_row); % number of row elements in merged subplot   

merged_height = subplot_rows*( height + gap_h ) - gap_h;   % merged subplot height
merged_width  = subplot_cols*( width  + gap_w ) - gap_w;   % merged subplot width

merged_bottom = (m-max(subplot_row)) * (height + gap_h) + marg_h(1);          % merged subplot bottom position
merged_left   = (min(subplot_col)-1) * (width  + gap_w) - gap_w + marg_w(1);  % merged subplot left position
pos_vec = [merged_left merged_bottom merged_width merged_height];

h_subplotv=subplot('Position', pos_vec, varg{:} );

if nargout > 0
   h_subplot = h_subplotv;
end

end

