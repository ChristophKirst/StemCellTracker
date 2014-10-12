function frmt = imfrmtFormatFromCell(data)
%
% format = imfrmtFormatFromCell(data)
%
% description:
%     determines the format of the cell array
% 
% input:
%     data    cell array
%
% output:
%     frmt    the cell format of the data
%
% note: 
%     take subset of 'SUVW' or '' if larger dims


d = ndims1(data);

if d <= 4
   frmt = 'SUVW';
   frmt = frmt(1:d);
else
   frmt = '';
end

end
      
   




