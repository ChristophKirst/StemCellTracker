function t = birc_getAquisitionTime(filename, format)
%
% t = bric_getAquisitionTime(filename, format)
%
% description: 
%    returns acquisition-time-local property of image taken by Bunty microscope in BRIC
%
% input:
%    filename   name of image file or folder
%    format     (optional) 'str'='string' or 'num', 'vec', 'seconds' ('string')
%
% output:
%    t          time as string or cell of strings
%
% See also: imfinfo

if nargin < 2
   format = 'string';
end

t = mmaquisitiontime(filename, format);

end
