function m = numpyToMat(np)
%
% m = numpyToMat(np)
%
% description:
%    converts numpy array to matalb array
%
% input:
%    np   numpy array
%
% output:
%   m     matlab array
%
% note: 
% nonly works for double so far, could extend to different data typesdec2

% m = cell2mat(cell(np.flatten.tolist));
m = typecast(uint32(char(np.tobytes)), 'single')

m = uint8(char(np.tobytes));

size(m)



sh = cell2mat(cell(np.shape));
m = reshape(m, sh);

end