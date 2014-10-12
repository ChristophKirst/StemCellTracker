function img = imfrmtPermute(img, infrmt, outfrmt)
%
% img = imfrmtPermute(img, in_format, out_format)
%
% description: 
%     takes a image and reformats it from informat into outformat
%     checks if no dimensions are lost or added
%
% input:
%     img        input image
%     infrmt     (optional) format of input image (imfrmtFormat(img))
%     outfrmt    (optional) format of output image 
%
% output:
%     img       reformatted image
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


if length(unique(lower(infrmt))) ~= length(infrmt)
   error('imfrmtReformat: labels appear more than once in input format: %s', infrmt)
end
if length(unique(lower(outfrmt))) ~= length(outfrmt)
   error('imfrmtReformat: labels appear more than once in output format: %s', outfrmt)
end

outfrmtl = lower(outfrmt);
infrmtl = lower(infrmt);

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
   if isempty(sf) % remove dimensions only if it is of size 1 (squeeze)
      if size(img,i) == 1
         outfrmtl = [outfrmtl, infrmtl(i)]; %#ok<AGROW> %% appending it effectively removes dimensions of size 1
      else
         error('imfrmtPermute: non-trivial dimension %s not specified in output format!', infrmtl(i), outfrmtl);
      end
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
%per

% premute
img = permute(img, per);

end
