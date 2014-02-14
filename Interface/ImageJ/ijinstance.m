function mij = ijinstance()
%
% mij = ijinstance()
%
% description: 
%     tries to find running ImageJ instance returns [] if not found
%
% output:
%     mij      ImageJ object
%
% See also: ijstart

try
   mij = ij.IJ.getInstance();
catch
   mij = []; 
end

end
   


