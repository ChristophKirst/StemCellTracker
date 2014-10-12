function img = imfrmtReformat(img, infrmt, outfrmt)
%
% img = imfrmtReformat(img, in_format, out_format)
%
% description: 
%     takes a image and reformats it from informat into outformat
%     removes all extra dimensions in img by using first index
%     adds singlentons in missing dimensions
%
% input:
%     img        input image
%     infrmt     (optional) format of input image (imfrmtFormat(img))
%     outfrmt    (optional) format of output image 
%
% output:
%     img        reformatted image
%
% note:
%     upper case letters are mathematically in postive axis, lower case mathematically inverted axis

if nargin < 2 || isempty(infrmt)
   infrmt = imfrmtFormat(img);
end
infrmt = imfrmtFormat(infrmt);

if nargin < 3 || isempty(outfrmt)
   outfrmt = 'XYZCT';
end
outfrmt = imfrmtFormat(outfrmt);

if strcmp(infrmt, outfrmt)
   return
end

if length(unique(lower(infrmt))) ~= length(infrmt)
   error('imfrmtReformat: labels appear more than once in input format: %s', infrmt)
end
if length(unique(lower(outfrmt))) ~= length(outfrmt)
   error('imfrmtReformat: labels appear more than once in output format: %s', outfrmt)
end

outfrmtl = lower(outfrmt);
infrmtl = lower(infrmt);
nout = length(outfrmt);

%reverse coordinates
[id, pos] = ismember(infrmtl, outfrmtl);
pos = pos(id);
id = find(id);

for i= 1:length(id)
   if infrmt(id(i)) ~= outfrmt(pos(i))
      img = flip(img, id(i));
   end
end

% add extra dimensions
for i = 1:length(outfrmtl)
   sf = strfind(infrmtl, outfrmtl(i));
   if isempty(sf)
      infrmtl(end+1) = outfrmtl(i);  %#ok<AGROW> % add to the end as this would introduce a new dimensions of size 1 in the original array
   end
end


% find dimensions to squeeze 
for i = 1:length(infrmtl)
   sf = strfind(outfrmtl, infrmtl(i));
   if isempty(sf) % remove dimensions
      outfrmtl = [outfrmtl, infrmtl(i)]; %#ok<AGROW> %% appending it effectively removes dimensions of size 1
   end
end
% now out_format and in_format have same size
%infrmt
%outfrmt

%find permutation
per = ones(1, length(infrmtl));
for i = 1:length(infrmtl)
   per(i) = strfind(infrmtl, outfrmtl(i));
end

if length(per) == 1
   per = [per,2];
end

% premute
img = permute(img, per);

%cut out reduced dimensions
asgn = repmat({':'}, 1, length(outfrmtl));
asgn(nout+1:length(outfrmtl)) = num2cell(ones(1,length(outfrmtl)-nout));
img = img(asgn{:});

end