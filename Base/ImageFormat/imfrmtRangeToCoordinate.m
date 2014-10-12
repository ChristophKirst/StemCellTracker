function coords = imfrmtRangeToCoordinate(isize, ifrmt, rgs)
%
% coords = imfrmtRangeToCoordinate(isize, ifrmt, trg)
%
% description: 
%     converts range to array of coordinates
%
% input:
%     isize      the size of the data array
%     ifrmt      the reference format
%     rgs        coordinate ranges
%
% output:
%     coords     coordinates
% 
% See also: imfrmtIndexToCoordinate

% generate indices
idx = imfrmtRangeToIndex(isize, ifrmt, rgs);
coords = imfrmtIndexToCoordinate(isize, ifrmt, idx);

end








   