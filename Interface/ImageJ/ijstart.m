function mij = ijstart()
%
% ijstart()
%
% description:
%    starts ImageJ program
%
% output:
%    mij  ImageJ instance
%
% See also: ijinitialize

mij = ijinstance();

if isempty(mij)
   try
      mij = ij.ImageJ([], 2);
   catch
      error('ijstart: could not start ImageJ, run ijinitialize first');
   end
end


