function [im1, bckgnd] = subtractBackground(im0)








end




% 
% %% Aryehs code 
% % imclose seems not the best choice
% nucAreaBckgnd = 100;
% 
% % filter a bit first to remove extreme local min, which otherwise will
% % cause holes around large bright objects after imclose
% hg = fspecial('gaussian', 12, 2);
% bckgnd = imfilter(im0, hg, 'replicate');
% 
% rdisk = ceil(sqrt(userParam.nucAreaBckgnd/pi));
% % remove the nuclei from image before filtering.
% bckgnd = imclose(bckgnd, strel('square', 2*rdisk+1));
% hg = fspecial('gaussian', 12*rdisk, 2*rdisk);
% bckgnd = imfilter(bckgnd, hg, 'replicate');
% mm = min(bckgnd(:));
% bckgnd = bckgnd - mm;
% im1 = im0 - bckgnd;

