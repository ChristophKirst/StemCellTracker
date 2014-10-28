function cols = detectColoniesFromAlignment(algn, varargin)
%
% cols = detectColoniesFromAlignment(algn, param)
%
% descritption:
%     detects shapes in the resamples aligned preview image
%     rescales rois to full scale and creates an array of colony classes
%
% input:
%    img    image
%    scale  resmapling scale
%    param  parameter struct
%           .threshold    threshold ([] = thresholdFirstMin(img))
%           .filtersize   size of Gaussian filter ([] = no filter) 
%           .strel        structure element for morphological closing ([] = strel('disk', 20))
%           .output       'ROIs' or 'label' ('label')
%           .plot         plot result
%           options for points2shape
%           
% output:
%    img    labeled image or array of ROIPolygon objects



% generate preview image



end