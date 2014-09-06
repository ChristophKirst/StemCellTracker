function bad = imcheckimage(img, warn, err)
%
% bad = imcheckimages(img, warn, err)
%
% description:
%     checks image for zero 2d subimages in image
%
% input:
%     img    image to check
%     warn   (optional) true to print warings ([] = false)
%     err    (optional) true to create error ([] = false)
% 
% output:
%     bad    array of bad images

% check for empty images

if nargin < 2
   warn = true;
end
if nargin < 3
   err = false;
end

format = imformat(img);

if isempty(format)
   bad = true;
   msg = 'imcheckimages: not valid image format';
   if warn, warning(msg); end
   if err, error(msg); end
end

zct = format(3:end);

switch length(zct)
   case 0
      bad = 0;
      if max(max(img)) == 0
         bad = true;
         msg = 'imcheckimages: zeros in 2d image';
         if warn, warning(msg); end
         if err, error(msg); end
      end
   case 1
      bad = zeros(size(img,3), 1);
      for z = 1:size(img,3)
         if max(max(img(:,:,z))) == 0
            bad(z) = true;
            msg = sprintf('imcheckimages: zero image at %s=%g', zct(1), z);
            if warn, warning(msg); end %#ok<SPWRN>
            if err, error(msg); end %#ok<SPERR>
         end
      end
   case 2
      bad = zeros(size(img,3), size(img,4));
      for z = 1:size(img,3)
         for c = 1:size(img,4)

            if max(max(img(:,:,z,c))) == 0
               bad(z,c) = true;
               msg = sprintf('imcheckimages: zero image at (%s,%s)=(%g,%g)', zct(1), zct(2), z, c);
               if warn, warning(msg); end %#ok<SPWRN>
               if err, error(msg); end %#ok<SPERR>
            end
         end
      end
   case 3
      bad = zeros(size(img,3), size(img,4));
      for z = 1:size(img,3)
         for c = 1:size(img,4)
            for t = 1:size(img,5)
               if max(max(img(:,:,z,c,t))) == 0
                  bad(z,c,t) = true;
                  msg = sprintf('imcheckimages: zero image at (%s,%s,%s)=(%g,%g,%g)', zct(1), zct(2), zct(3), z, c, t);
                  if warn, warning(msg); end %#ok<SPWRN>
                  if err, error(msg); end %#ok<SPERR>
               end
            end
         end     
      end
end

end