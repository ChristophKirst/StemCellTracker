function rgs = imfrmtReformatRange(inSize, inFrmt, outFrmt, rgs)
%
% img = imfrmtReformatRange(inSize, inFrmt, outFrmt, tgr)
%
% description: 
%     takes coordinate ranges rgs and reformats its entries to outfrmt assuming image of size isize
%     extra dimensions in outfrmt are not added
%     all entries in rgs not specified in outfrmt are removed
%
% input:
%     inSize     size of input data        
%     inFrmt     format of input data
%     outFrmt    the output format
%     rgs        coordinate range of the input data 
%
% output:
%     rgs        reformatted ranges
%
% note: 
%     does not reformat order of coord names in the struct
% 
% See also: imfrmtAssignment

inFrmtl = lower(inFrmt);

tgrNames = fieldnames(rgs);
tgrNamesl = lower(tgrNames);
outFrmtl = lower(outFrmt);

%find all entries not in output format and delete
delfrmt = setdiff(inFrmtl, outFrmtl);

id = ismember(tgrNamesl, num2cell(delfrmt));
rgs = rmfield(rgs, tgrNames(id));

% find entries to invert
tgrNamesl = tgrNamesl(~id);
tgrNames = tgrNames(~id);

[id, pos] = ismember(tgrNamesl, num2cell(outFrmtl));
id = find(id);

% invert if requested
for i = 1:length(id)
   if tgrNames{id(i)} ~= outFrmt(pos(i))
      o = outFrmtl(pos(i));
      k = find(inFrmtl == o, 1, 'first');
      o = outFrmt(pos(i));
      rgs = renameStructField(rgs, tgrNames{id(i)}, o);
      rgs.(o) = inSize(k) - rgs.(o) + 1;
   end
end


end
   