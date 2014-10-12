function asgn = imfrmtAssignment(inSize, inFrmt, varargin)
%
% asgn = imfrmtAssignment(isize, ifrmt, dataSpecs)
%
% description:
%     returns cell asgn such that for an image of size isize and format ifrmt
%     img(asgn{:}) is given by the data specificatons dataSpecs
%     
% 

param = parseParameter(varargin);
pnames = fieldnames(param);

d = length(inFrmt);
asgn = repmat({':'}, 1, d);

[id, pos] = ismember(lower(pnames), num2cell(lower(inFrmt)));
pos = pos(id);
id = find(id);

for i = 1:length(id)
   vals = param.(pnames{id(i)});
   if iscell(vals)
      vals = cell2mat(vals);
   end
   if pnames{id(i)} == inFrmt(pos(i))
      asgn{pos(i)} = vals;
   else
      asgn{pos(i)} = inSize(pos(i)) - vals + 1;
   end
end

end