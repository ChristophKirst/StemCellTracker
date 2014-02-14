function id = imarisid(imaris)
%
% id = imarisid(imaris)
%
% description:
%    returns imaris application id
%
% input:
%    imaris   (optional) Imaris application reference


id = [];

if nargin < 1 
   try 
      ilib = ImarisLib();
      server = ilib.GetServer();
      id = server.GetObjectID(0);
   catch
   end
   
elseif isimaris(imaris)
   % TODO
   
end
