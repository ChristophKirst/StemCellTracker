function coords = imfrmtIndexToCoordinate(isize, ifrmt, id)
%
% coords = imfrmtIndexToCoordinate(isize, ifrmt, id)
%
% description: 
%     converts index to an array of coordinates struct given the size and format and the indices
%
% input:
%     isize      the size of the data array
%     ifrmt      the reference format
%     id         the coordinate index
%
% output:
%     coords     (array of) coordinates
% 
% See also: imfrmtReformatRange


idx = imind2sub(isize, id);
idx = num2cell(idx);

nfrmt = length(ifrmt);
args = cell(1, 2*nfrmt);
for i = 1:nfrmt
   args{2*i-1} = ifrmt(i);
   args{2*i}   = idx(:,i)';
end

coords = struct(args{:});

end
   