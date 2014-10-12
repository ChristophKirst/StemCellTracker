function frmt = imfrmtRemoveFormat(frmt, rmFrmt)
%
% frmt = imfrmtRemoveFormat(frmt, rmFrmt)
%
% description: 
%     returns data format obtained by removing rmFrmt 
%
% input:
%     frmt        format of data
%     rmFrmt      formats to remove
%
% output:
%     frmt        data format with removed entries
%
% note:
%     lower/ upper case invariant


frmtl = lower(frmt);
rmFrmtl = lower(rmFrmt);

id = ismember(frmtl, rmFrmtl);
frmt = frmt(~id);

end