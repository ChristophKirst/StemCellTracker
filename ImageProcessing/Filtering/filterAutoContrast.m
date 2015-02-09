function img = filterAutoContrast(img)
%
% img = filterAutoContrast(img)
%
% description:
%   Auto contrast image such that ideally objects
%   have intensity range ]0.5 1] and background  [0 0.5]
%
% input:
%   img      image, assumed to be in rage 0 1
%
% output:    
%   img      img after auto-contrast.
%

% Based on: autocontrast by Aaron Ponti, 2005/08/30

if max(img(:)) > 1
   warning('filterAutoContrast: image not in range 0..1, clipping the image!');
   img = imclip(img, 0, 1);
end

% detect limits for imadjust
minTol=0.02;
lim=[0 1];
while lim(1)==0 && lim(2)==1
    minTol=2*minTol;
    maxTol=1-minTol;
    if minTol>=maxTol
        break;
    end
    lim = stretchlim(img, [minTol maxTol]);
end

img = imadjust(img,lim);

end