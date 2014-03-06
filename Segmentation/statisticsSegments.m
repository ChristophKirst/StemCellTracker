function stats = statisticsSegments(image, label, props)
%
% stats = statisticsSegments(image, label, props)
%
% description:
%      generates statistics for labels in labeled image label using the underlying image image
%
% input:
%      image   image to measure from
%      label   labeled image indicating regions to calc statistics for
%      props   (optional) cell array of names of statistical properties to calculate ([] = {} = 'all')
%
% output:
%      stats   struct array with results
%
% See also: regionprops

