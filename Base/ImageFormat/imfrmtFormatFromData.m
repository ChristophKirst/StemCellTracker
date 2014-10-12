function frmt = imfrmtFormatFromData(data)
%
% format = imfrmtFormatFromData(data)
%
% description:
%     determines the format of the data
%     channel dimension is detected as dimensions having size <= 3
% 
% input:
%     data    image data
%
% output:
%     frmt    the data format of the data
%
% note: 
%     inferring data form an array is ambiguous:
%     'XYZ' = 'XYT', 'XYCZ' = 'XYCT', 'XYZC' = 'XYTC'

frmt = imfrmtFormatFromSize(size(data));

end
      
   




