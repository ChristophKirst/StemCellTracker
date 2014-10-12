function [si, frmt] = imreadBFSize(name, frmt)
%
% pos = imreadBFSize(name)
%
% decription:
%    returns sizes for bio-format image

% input:
%    name    filename or reader
%    frmt    (optional) format ('XYZCTS')
%             
% output:
%    si      size in XYZCTS format 
%    frmt    format specified in input

ireader = imreadBFReader(name);
x = ireader.getSizeX();
y = ireader.getSizeY();
z = ireader.getSizeZ();
c = ireader.getSizeC();
t = ireader.getSizeT();
s = ireader.getSeriesCount();

si = [x,y,z,c,t,s];

if nargin < 2
   frmt = 'XYZCTS';
end

si = imfrmtReformatSize(si, 'XYZCTS', frmt);

end