function mij = ijinstance(startij)
%
% mij = ijinstance()
%
% description: 
%     tries to find running ImageJ instance
%
% output:
%     mij      ImageJ object
%     startij  start imagej if no instance is found (0)
%
% See also: ijstart

if nargin < 1
   startij = 1;
end

mij = ij.IJ.getInstance();

if isempty(mij) && startij
   mij = ijstart();
end
   


