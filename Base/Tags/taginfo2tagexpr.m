function texpr = taginfo2tagexpr(tsplit, tinfo)
%
% texpr = taginfo2tagexpr(tsplit, tinfo)
%
% description:
%    generates a tagexpr string from the splitting and tifno struc obtained by tagexpr2tagnames
%
% input:
%    tplsit      split of the tagexpr
%    tinfo       struct with the tag info
%
% output:
%    tfrmt       the tag format 
%
% See also: tagexpr2tagnames

tnames = {tinfo.name};

res = repmat({''}, length(tsplit)-1, 1);
for i = 1:length(tnames)
   res(tinfo(i).pos) = tinfo(i).tag;
end

texpr = tsplit{1};
for i = 1:length(res);
   texpr = [texpr,  res{i}, tsplit{i+1}]; %#ok<AGROW>
end

end
   
   



