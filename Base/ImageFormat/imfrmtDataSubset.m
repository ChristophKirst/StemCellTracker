function data = imfrmtDataSubset(data, inFrmt, varargin)
%
% data = imfrmtDataSubset(data, ifrtm, varargin)
%
% description:
%    extract a subset of the data given the data specifications datasepc
%
% See also: imfrmtAssignment

asgn = imfrmtAssignment(size(data), inFrmt, varargin);
data = data(asgn{:});
 
end