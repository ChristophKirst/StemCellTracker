function np = numpyFromMat(m)
%
% m = numpyToMat(np)
%
% description:
%    converts numerical array m to numpy array np
%
% input:
%    m    matlab array
%
% output:
%    np   numpy array

np = py.numpy.array(m(:)');

si = size(m);
si(1:2) = si(2:-1:1);
np = np.reshape(si);

end