function rgb3d = label2rgb3d(varargin)
%
% label2rgb3d(label, camp, zerocolor, order)
%
% description:
%     converts labeled 3d image into colorized 3d image using a color map
%     works as label2rgb for 3d images
%
% input: 
%     label     labeled 3d image
%     cmap      (optional) color map as in label2rgb
%     zerocolor (optional) color for background
%     order     (optional) 'shuffle' or 'noshuffle'
% 
% output:
%     rgb3d      3d color image 
%
% notes: based on label2rgb
%
% See also: label2rgb, imcolorize


[label,map,zerocolor,order,fcnflag] = parse_inputs(varargin{:});

% Determine the number of regions in the label matrix.
numregion = double(max(label(:)));

% If MAP is a function, evaluate it.  Make sure that the evaluated function
% returns a valid colormap.
if  fcnflag == 1
    if numregion == 0
      cmap = [];
    else
      cmap = feval(map, numregion);
      if ~isreal(cmap) || any(cmap(:) > 1) || any(cmap(:) < 0) || ...
            ~isequal(size(cmap,2),3) || size(cmap,1) < 1
        error(message('images:label2rgb:functionReturnsInvalidColormap'));
      end
    end
else
    cmap = map;
end

% If ORDER is set to 'shuffle', create a private stream with a fixed seed,
% which creates the same "random" permutation every time it is called.
if isequal(order,'shuffle')
    stream = RandStream('swb2712','seed',0);
    index = randperm(stream,numregion);
    cmap = cmap(index,:,:);
end

% Issue a warning if the zerocolor (boundary color) matches the color of one
% of the regions. 
for i=1:numregion
  if isequal(zerocolor,cmap(i,:))
    warning(message('images:label2rgb:zerocolorSameAsRegionColor', i));
  end
end
cmap = [zerocolor;cmap];


% Make sure A is in the range from 1 to size(cm,1)
label = max(1,min(label,size(cmap,1)));

% Extract r,g,b components
r = zeros(size(label)); r(:) = cmap(label,1);
g = zeros(size(label)); g(:) = cmap(label,2);
b = zeros(size(label)); b(:) = cmap(label,3);

if nargout==3,
  rgb3d = r;
else
  rgb3d = zeros([size(r),3]);
  rgb3d(:,:,:,1) = r;
  rgb3d(:,:,:,2) = g;
  rgb3d(:,:,:,3) = b;
end

end




%  Function: parse_inputs
%  ----------------------
function [L, Map, Zerocolor, Order, Fcnflag] = parse_inputs(varargin) 
% L         label matrix: matrix containing non-negative values.  
% Map       colormap: name of standard colormap, user-defined map, function
%           handle.
% Zerocolor RGB triple or Colorspec
% Order     keyword if specified: 'shuffle' or 'noshuffle'
% Fcnflag   flag to indicating that Map is a function


%narginchk(1,4);

% set defaults
L = varargin{1};
Map = 'jet';    
Zerocolor = [1 1 1]; 
Order = 'noshuffle';
Fcnflag = 0;

% parse inputs
if nargin > 1
    Map = varargin{2};
end
if nargin > 2
    Zerocolor = varargin{3};
end
if nargin > 3
    Order = varargin{4};
end

% error checking for L
validateattributes(L,{'numeric','logical'}, ...
              {'real' 'nonsparse' 'finite' 'nonnegative' 'integer'}, ...
              mfilename,'L',1);

% error checking for Map
[fcn, fcnchk_msg] = fcnchk(Map);
if isempty(fcnchk_msg)
    Map = fcn;
    Fcnflag = 1;
else
    if isnumeric(Map)
        if ~isreal(Map) || any(Map(:) > 1) || any(Map(:) < 0) || ...
                    ~isequal(size(Map,2), 3) || size(Map,1) < 1
          error(message('images:label2rgb:invalidColormap'));
        end
    else
        error(fcnchk_msg);
    end
end    
    
% error checking for Zerocolor
if ~ischar(Zerocolor)
    % check if Zerocolor is a RGB triple
    if ~isreal(Zerocolor) || ~isequal(size(Zerocolor),[1 3]) || ...
                any(Zerocolor> 1) || any(Zerocolor < 0)
      error(message('images:label2rgb:invalidZerocolor'));
    end
else    
    %[cspec, msg] = cspecchk(Zerocolor);
    msg = '';
    if ~isempty(msg)
	%message is translated at source.
        error(message('images:label2rgb:notInColorspec', msg))
    else
        %Zerocolor = cspec;
        Zerocolor = [0 0 0];
    end
end

% error checking for Order
valid_order = {'shuffle', 'noshuffle'};
idx = strncmpi(Order, valid_order,length(Order));
if ~any(idx)
    error(message('images:label2rgb:invalidEntryForOrder'))
elseif nnz(idx) > 1
    error(message('images:label2rgb:ambiguousEntryForOrder', Order))
else
    Order = valid_order{idx};
end

end

  