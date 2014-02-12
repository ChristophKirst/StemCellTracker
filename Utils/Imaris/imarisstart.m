function imaris = imarisstart()
%
% imaris = istart(ipath)
%
% description:
%    starts Imaris searching for XT in ipath
%
% input:
%    ipath   path to ImarisXT or id of instance TODO !
%
% 
% TODO: make this independent of IceImarisConnector
% See also: IceImarisConnector


imaris = imarisinstance();

if isempty(imaris)

    conn  = IceImarisConnector([],1);
    success = conn.startImaris();

    if ~success
        error('imarisstart: could not start Imaris');
    end

    imaris = conn.mImarisApplication;

end




