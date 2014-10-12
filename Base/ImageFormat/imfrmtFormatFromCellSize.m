function frmt = imfrmtFormatFromCellSize(csize)
%
% frmt = imfrmtFormatFromCellSize(csize)
%
% description:
%     determines the format of the cell array of size csize
% 
% input:
%     csize   size of cell array
%
% output:
%     frmt    the cell format of the data
%
% note: 
%     take subset of 'UVW' or '' if larger dims


d = length(csize);

if d <= 3
   frmt = 'SUVW';
   frmt = frmt(1:d);
else
   frmt = '';
end

end
      
   




