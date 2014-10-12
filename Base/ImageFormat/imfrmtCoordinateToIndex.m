function idx = imfrmtCoordinateToIndex(isize, ifrmt, coord)
%
% idx = imfrmtCoordinateToIndex(isize, ifrmt, coord)
%
% description: 
%     converts (array of) coordinates to index given the size and format
%
% input:
%     isize      the size of the data array
%     ifrmt      the reference format
%     coord      coordinate or array of coordinates
%
% output:
%     id         array of indices
% 
% See also: imfrmtIndexToCoordinate


% generate indices

nfrmt = length(ifrmt);
nids  = length(coord);
sub = zeros(nids, nfrmt);

for i = 1:nfrmt
   sub(:,i) = [coord.(ifrmt(i))];
end

idx = imsub2ind(isize, sub);

end
   