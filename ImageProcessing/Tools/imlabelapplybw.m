function labela = imlabelapplybw(label, fun)
%
% labela = imlabelapplybw(label, fun)
%
% description:
%      applies fun that transfroms a bw image into a bw onto each label separately
%
% input:
%      label     labeled image
%      fun       function operating on bw image representing a single region
% 
% output:
%      labela    results of applying fun
%
% See also: imlabel

lab = imlabel(label);
labela = label;

for l = lab
   imgl = fun(label == l);
   labela(imgl > 0) = l;
end

end