function frmt = imfrmtFormat(in)
%
% frmt = imfrmtFormat(frmt)
%
% description:
%     determines image format of the input in 
%
% input:
%     in    data array, size array (size(in,1) == 1), cell array, format, ImageSource, ImageInfo
%
% output:
%     frm   standarized format

if isnumeric(in) || islogical(in)
   if size(in,1) == 1 && ismatrix(in);
      frmt = infrmFormatFromSize(in);
   else
      frmt = imfrmtFormatFromData(in);
   end
   
elseif iscell(in)
   frmt = imfrmtFormatFromCell(in);
   
elseif ischar(in)
   frmt = in;
   
   switch frmt
      case 'matlab'
         frmt = 'yXCZT';  
      case 'imagej'
         frmt = 'XYCZT';  % to check
      case 'data'
         frmt = 'XYZCT';
      case 'cell'
         frmt = 'SUVW';
   end
   
elseif isa(in, 'ImageSource') || isa(in, 'ImageInfo')
   frmt = in.dataFormat;
   
else
   error('imfrmtFormat: not a valid input !');
end

end

