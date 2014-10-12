function id = imfrmtReformatIndex(id, isize, infrmt, outfrmt)
%
% id = imfrmtReformatSize(id, isize infrmt, outfrmt)
%
% description: 
%     takes indices of an image of size isize anf format informt and
%     reformats them to outfrmt
%
% input:
%     isize      input image size
%     infrmt     (optional) format of input image (imfrmtFormat(img))
%     outfrmt    (optional) format of output image 
% output:
%     id         reformatted image size
%

sub = imind2sub(isize, id);
[sub, isize] = imfrmtReformatSubIndex(sub, isize, infrmt, outfrmt);
id = imsub2ind(isize, sub);

end
   