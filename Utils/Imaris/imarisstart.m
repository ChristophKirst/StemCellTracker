function imaris = imarisstart(ipath)
%
% imaris = istart(ipath)
%
% description:
%    starts Imaris searching for XT in ipath
%
% input:
%    ipath   path to ImarisXT or id of instance
%
% 
% TODO: make this independent of IceImarisConnector
% See also: IceImarisConnector

conn  = IceMarisConnector([],1);
success = conn.startImaris();

if ~success
   error('imarisstart: could not start Imaris');
end

imaris = conn.mImarisApplication;

end




