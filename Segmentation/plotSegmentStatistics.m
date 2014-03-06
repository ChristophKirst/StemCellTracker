function stats = plotSegmentStatistics(image, label, statnames)
%
% stats = plotSegmentStatistics(stats)
% stats = plotSegmentStatistics(image, label, statnames)
%
% description:
%     plots segment statistics of the labeled image label for underlying image
%
% input:
%     image     underlying image to measure from
%     label     labeled image
%     statnames 
%     stats     struct array with segment properties
%
% output:
%     

if nargin == 1
   stats = image;
elseif naragin == 2
   stats = statisticsSegments(image, label);
else
   stats = statisticsSegments(image, label, statnames);
end

statnames = fieldnames(stats);
% plot a histogram and a cumulative distribution for the proeprties



