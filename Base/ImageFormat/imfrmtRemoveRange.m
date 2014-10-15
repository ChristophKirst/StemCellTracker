function range = imfrmtRemoveRange(range, rmFrmt)
%
% frmt = imfrmtRemoveFormat(frmt, rmFrmt)
%
% description: 
%     returns range obtained by removing entries rmFrmt (case insensitive)
%
% input:
%     range       data range
%     rmFrmt      formats to remove
%
% output:
%     frmt        data format with removed entries
%
% note:
%     lower/ upper case invariant

frmt = fieldnames(range);
frmtl = lower(frmt);
rmFrmtl = lower(rmFrmt);
if ischar(rmFrmtl)
   rmFrmtl = num2cell(rmFrmtl);
end

id = ismember(frmtl, rmFrmtl);
range = rmfield(range, frmt(id)) ;

end