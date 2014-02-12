function imaris = imarisinstance(imarisApplicationID)
%
% imaris = imarisinstance(imarisApplicationID)
%
% description:
%    returns a running Imaris application instance
%
% input:
%    imarisApplicationID    (optional) id of the Imaris application
%
% See also: imarisstart

try
   if nargin < 1
      vImarisLib = ImarisLib();
      server = vImarisLib.GetServer();
      imarisApplicationID = server.GetObjectID(0);
      imaris = vImarisLib.GetApplication(imarisApplicationID);
   elseif isa(imarisApplicationID, 'Imaris.IApplicationPrxHelper')
      imaris = imarisApplicationID;
   else
      if ischar(imarisApplicationID)
         imarisApplicationID = round(str2double(imarisApplicationID));   
      end
      vImarisLib = ImarisLib();
      imaris = vImarisLib.GetApplication(imarisApplicationID);
   end
   
catch %#ok<CTCH>
   %error('imarisintance: no instance running, path to ImarisLib not set, or wrong id');
   imaris = [];
end
