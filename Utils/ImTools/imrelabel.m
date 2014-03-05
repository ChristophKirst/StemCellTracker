function relabel = imrelabel(label)
%
% relabel = imrelabel(label)
%
% description:
%    relabels the labeled image using labels 1 to number of labeled regoins
%
% input:
%    label    labeled image (2D/3D)
%
% output:
%    relabel  relabeld image
%
% See also: imlabelseparate

labs = imlabel(label);
relabel = label;

lnew = 1;
for lold = labs
   relabel(label == lold) = lnew;
   lnew = lnew + 1;
end

end
   