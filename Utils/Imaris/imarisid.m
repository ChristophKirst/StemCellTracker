function id = imarisid()
%
% id = imarisid()
%
% description:
%    returns imaris application id

vImarisLib = ImarisLib();
server = vImarisLib.GetServer();
id = server.GetObjectID(0);

end
