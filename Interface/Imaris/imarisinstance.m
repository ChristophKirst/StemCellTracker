function imaris = imarisinstance(id)
%
% imaris = imarisinstance(id)
%
% description:
%    returns a running Imaris application instance
%    or [] if a valid Imaris application cannot be found
%
% input:
%    id    (optional) id of a specific Imaris application
%
% See also: imarisstart

imaris = [];

try
   ilib = ImarisLib();
   if isempty(ilib)
      return
   end
   
   server = ilib.GetServer();
   if isempty(server)
      return
   end
   
   napps = server.GetNumberOfObjects();
   if napps == 0
      return
   end
    
   if nargin < 1
      id = server.GetObjectID(0);
   elseif ischar(id)
      id = round(str2double(id));
   end
   
   for a = 0 : napps - 1
      if server.GetObjectID(a) == id
         imaris = ilib.GetApplication(id);
         return
      end
   end
catch 
end

end
   
