function imglab = polygonToLabeledImage(pol, varargin)
%
% imglab = polygonToLabeledImage(pol, varargin)
%
% description:
%     converts a polygon into a labeled image
%
% input:
%     pol     polygon
%     param   parameter struct with entries
%             size     image size ([] = bounding box of the poly)
%             range    coords of the image as [lower, upper] ([] = bounding box of the poly)
%             
% output:
%     imglab  labled image
%
% See also: polygonFromLabeledImage, polygonFromMask

param = parseParameter(varargin);

si = getParameter(param, 'size', []);
rg = getParameter(param, 'range', []);


if isempty(rg)
   bb = polygonToBoundingBox(pol);
   rg = sort(bb(:, [1,3]),2);
end

if isempty(si)
   si = diff(rg');
   si(si < 1) = 1;
end

% scale the polygons appropiately    
   
   


BW = false([m,n]);

% Start with any enclosing boundaries not enclosed by others
[rr,~] = find(A);
listToAdd = false(1,length(XY));
listToAdd(setdiff(1:length(XY),rr)) = true;
while any(listToAdd)
    % Add the next enclosing mask
    nextAdd = find(listToAdd,1);
    BW = BW | pMaskFcn(XY{nextAdd});
    listToAdd(nextAdd) = false;
    % Subtract any enclosed masks (and then add their children)
    for nextSub = find(A(:,nextAdd))'
        BW = BW & ~pMaskFcn(XY{nextSub});
        % Add any children 
        listToAdd(A(:,nextSub)) = true;
    end
end

function [XY, A, pMaskFcn, m, n] = parseInputs(XY, varargin)

% Ensure cell input
if isnumeric(XY)
    nanIdxs = [0; find(any(isnan(XY),2)); size(XY,1)];
    XYcell = cell(length(nanIdxs)-1,1);
    for i = 1:length(XYcell)
        inds = nanIdxs(i)+1 : nanIdxs(i+1)-1;
        if isempty(inds)
            XYcell{i} = zeros(0,2);
        else
            XYcell{i} = XY(inds,:);
        end
    end
    XY = XYcell;
end
numXY = length(XY);

% Detect IJ style if offered in params/values
if nargin>2 ...
        && ischar(varargin{end})   && strcmpi(varargin{end},'ij') ...
        && ischar(varargin{end-1}) && strcmpi(varargin{end-1},'style')
    XYcols = [2 1];
    varargin(end-1:end) = [];
else
    XYcols = [1 2];
end

% Detect association matrix or default to all-union
if islogical(varargin{end})
    A = varargin{end};
    varargin(end) = [];
else
    A = sparse([],[],[],numXY,numXY,false);
end

% If we are left with two vectors, we have XVEC, YVEC input. If we are left
% with a size matrix, we have BWSIZE input. Otherwise error.
if length(varargin)==1 && length(varargin{1})==2
    m = varargin{1}(1);
    n = varargin{1}(2);
    pMaskFcn = @(xy)poly2mask(xy(:,XYcols(1)),xy(:,XYcols(2)), m, n);
elseif length(varargin)==2 && isvector(varargin{1}) && isvector(varargin{2})
    xVec = varargin{1};
    yVec = varargin{2};
    m = length(yVec);
    n = length(xVec);
    pMaskFcn = @(xy)poly2mask(...
        interp1(xVec, 1:n, xy(:,XYcols(1)),'linear','extrap'),...
        interp1(yVec, 1:m, xy(:,XYcols(2)),'linear','extrap'),...
        m, n);
else
    error('mpoly2mask:input','input parameters could not be resolved.')
end