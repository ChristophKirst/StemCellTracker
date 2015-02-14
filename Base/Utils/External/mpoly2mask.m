function BW = mpoly2mask(XY, varargin)
% MPOLY2MASK Convert multiple region-of-interest polygons to a mask.
%
% BW = mpoly2mask(XY, BWSIZE) computes a binary region-of-interest mask,
% BW, from multiple region-of-interest polygons represented by XY. The size
% of BW (in rows, columns) is given in the 2-element BWSIZE. The class of
% BW is logical with 1 inside the set of XY polygons and 0 outside.
%
% XY is a cell array with separate polygon (xy) coordinates given as an
% N-by-2 array in each cell element. Alternatively, XY can be an N-by-2
% array containing all polygon coordinates, with each successive polygon
% separated by a pair of NaN elements. By default, the output mask is the
% union (mask1 & mask2 & ...) of all separate polygons.
%
% BW = mpoly2mask(..., 'style','ij') will interpret the input given in XY
% as being IJ-style contours rather than XY-style. This is consistent with
% MATLAB's bwboundaries command which returns contours in ij format (row
% coordinates in the first column, column coordinates in the second).
%
% BW = mpoly2mask(XY, XVEC, YVEC) where XVEC and YVEC are vectors,
% specifies the locations of the pixel centers of BW
%
% BW = mpoly2mask(XY, ..., A) allows complex relationships between the
% contours in XY. A is a square logical matrix with side length equal to
% the number of contours in XY, whose rows and columns correspond to the
% position of those contours in XY. The boundaries enclosed by the (i)th
% contour as the boundary enclosing the i(th) contour can both be found
% using A as follows:
% 
%    enclosing_boundary  = find(A(i,:));
%    enclosed_boundaries = find(A(:,i));
%
% For example, a donut shape can be represented by two circular contours in
% XY (XY{1} being the larger contour, X{2} being the smaller, and A having
% contents:
%
%     A = [0 0
%          1 0]
%
% Here, A(2,1) indicates that the second circle in XY is enclosed by the
% first. A may be sparse, as created by the command sparse(2,1,true,2,2).
% The logical matrix A can be ommited, in which case it defaults to a fully
% false matrix indicating no enclosing contours, and the contents of BW
% will be the union of all contours in XY.
%
% Example 1:
%     BW = imread('blobs.png');
%     [B,~,~,A] = bwboundaries(BW);
%     BW2 = mpoly2mask(B,size(BW),A,'style','ij');
%     figure
%     subplot(1,2,1), imshow(BW), title('Original')
%     subplot(1,2,2), imshow(BW2), title('mpoly2mask recreation')
%
% Note that mask with contours extracted using bwboundaries that are then
% directly passed to poly2mask (as in the example above) may have 1-pixel
% border regions that do not match. This issue, including a simple
% work-around is discussed at:
% http://blogs.mathworks.com/steve/2014/03/27/comparing-the-geometries-of-bwboundaries-and-poly2mask/

%%
[XY, A, pMaskFcn, m, n] = parseInputs(XY, varargin{:});
BW = false([m,n]);

% Start with any enclosing boundaries not enclosed by others
[rr,~] = find(A);
listToAdd = false(1,lengtmpoly2mask.mh(XY));
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