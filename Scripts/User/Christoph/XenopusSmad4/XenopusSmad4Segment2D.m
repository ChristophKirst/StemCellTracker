function [imgseg, imgmask, imgmed] = XenopusSmad4Segment2D(imgmp, verbose, time)

%% denoise
   
imgmed = filterMedian(imgmp, 5);
imgmed = imgmed - imopen(imgmed, strel('disk', 30));
imgmask = imgmed > 20;
imgmedmask = immask(imgmed, imgmask);

if verbose
   figure(1); clf
   implottiling({imgmp, imgmed; 255* imgmask, imgmedmask})
end



%% Seeding

imgf = filterDisk(imgmp, 6);

%stackf = imclip(stackdisk - mean(stackdisk(:)));
%stackf = stackf .* cast(stackmask, class(stackf));
%imgf = imgmed;
%imgf = filterDisk(imgf, 6);
%imgf = log(mat2gray(imgf)+eps) + 5;
%imgf(imgf <0) = 0;
%imgf = mat2gray(imgf);

imgmax = imextendedmax(mat2gray(imgf), 0.025);
%stackmax = imregionalmax(stackf);
imgmax = immask(imgmax, imgmask);

%
%stackmax = imerode(stackmax, strel('disk',3));
%stackmax = imdilate(stackmax, strel('disk',3));


if verbose
%    figure(10)
%    clf
%    set(gcf, 'Name', ['PreSeeding: ' filenames]);
%    is = imgf - mean(imgf(:));
%    is(is < 0) = 0;
%    implot3d(is);
%    
%    figure(11)
%    clf
%    set(gcf, 'Name', ['Seeding: ' filenames]);
%    implot3d(mat2gray(imgmax));
   
   % overlay in imagej
   
   figure(15); clf
   set(gcf, 'Name', ['Time: ', num2str(time), ' Seeding']);
   implottiling({imgf; imoverlay(mat2gray(imgmp), imgmax)});
   
end


%% Join Seeds

% [imgjoin, pairs, joins] = joinSeedsByRays(bwlabel(imgmax), mat2gray(imgmp),mat2gray(imgradient(imgmed)),...
%       'threshold.min', 0.2, ...
%       'threshold.max', 0.8,...
%       'threshold.change', 0.8,...
%       'threshold.gradient', [],...
%       'cutoff.distance', 20, ...
%       'averaging.ksize', 1,...
%       'plot.profiles', false);
% 
% 
% if verbose
%    figure(3); clf
%    implot(imoverlaylabel(imgmed, imgjoin))
%    
%    figure(4); clf
% %    imsubplot(1,2,1)
% %    implot(imoverlay(img, imgseg));
% %    for i = 1:size(pairs, 1);
% %       col = hsv2rgb([i/size(pairs,1), 1, 1]);
% %       plotSeedPairs(imgseg, pairs(i, :), col)
% %    end
%    
%    %^imsubplot(1,2,2)
%    implot(imoverlaylabel(imgmp, imgmax, false));
%    for i = 1:size(joins,1)
%       %col = hsv2rgb([i/size(joins,1), 1, 1])
%       plotSeedPairs(imgseg, joins(i,:),  [1,0,0])
%    end
% 
% end  


%% Segment


% 3d water shed on image + gradient 
%mi = max(imgmed(:));
%mg = max(imggrad(:));
%imgws = imgmed - 0.25 * mi / mg * imggrad;
%imgws(imgws < 0) = 0;

imgws = imgmed;
imgwsmin = imimposemin(iminvert(imgws), imgmax);

imgseg = watershed(imgwsmin);
imgseg = immask(imgseg, imgmask);



%% remove small segments

imgseg = postProcessSegments(imgseg, 'volume.min', 40,  'volume.max', 800, 'relabel', true);


if verbose
   
   figure(15); clf
   set(gcf, 'Name', ['Time: ', num2str(time), ' Watershed']);
   imgovl = imoverlaylabel(imgmp, impixelsurface(imgseg), false);
   implot(imgovl);
   
end


end

