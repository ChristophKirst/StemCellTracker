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
%
% note:
%    arrays are transformed such that m(x,y,...) = np[x-1, y-1,...]
%    in particular size(m) = np.shape



si = size(m);
m = permute(m, ndims(m):-1:1);

np = py.numpy.array(m(:)');
np = np.reshape(si);


% via filewriting is much slower
% tmpfn = [tempname, '.mat'];
% save(tmpfn, 'm')
% 
% np = py.scipy.io.loadmat(tmpfn);
% np = np.get('m');


%si(1:2) = si(2:-1:1);
%np = py.numpy.array(m(:)', pyargs('order', 'C'));

end