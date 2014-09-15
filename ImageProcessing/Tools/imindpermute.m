function ids = imindpermute(si, per, i)
%
% ids = imindpermute(si, per, i)
%
% description:
%    returns indices of an arrya a that has be transformed via b = permute(reshape(a,si), per)
%    then b(ids) = a(i)

ids = imind2sub(si, i);
ids = ids(:, per);
ids = reshape(imsub2ind(si(per), ids), size(i));

end