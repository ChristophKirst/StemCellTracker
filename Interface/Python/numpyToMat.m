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
%   as internalconversion is slow we create a temp file and write from there
%   converting to byte char object and then to matlab char looses characters

tmpfn = [tempname, '.mat'];
st = struct('data', np);
py.scipy.io.savemat(tmpfn, st);

m = load(tmpfn);
m = m.data;

delete(tmpfn);


% misc code
%m = cell2mat(cell(np.flatten.tolist));
%m = typecast(uint64(char(np.tobytes)), 'double')
%m = typecast(uint32(char(np.tobytes)), 'single')
%m = uint8(char(np.tobytes));
%sh = cell2mat(cell(np.shape));
%m = reshape(m, sh);

end