function imglab = segmentBySLIC(img, varargin)
%
% imglab = segmentBySLIC(img, varargin)
%
% description:
%    use Simple Linear Iterative Clustering in 5d color pixel space to determine pixels
% 
% input:
%    img    image, will be converted to unit8
%    param  parameter struct with entries
%           .superpixel   number of super pixels (2500)
%           .compactness   compactness factor (10)
%
% output:
%    imglab labled image of super pixels
%
% References: 
%    Radhakrishna Achanta, Appu Shaji, Kevin Smith, Aurelien Lucchi, Pascal Fua, and Sabine Süsstrunk, 
%    SLIC Superpixels Compared to State-of-the-art Superpixel Methods, 
%    IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 34, num. 11, p. 2274 - 2282, May 2012.
%
%    Radhakrishna Achanta, Appu Shaji, Kevin Smith, Aurelien Lucchi, Pascal Fua, and Sabine Süsstrunk, 
%    SLIC Superpixels, EPFL Technical Report no. 149300, June 2010.
%
%    http://ivrg.epfl.ch/research/superpixels
%
% See also: slicmex

param = parseParameter(varargin);

nspxls = getParameter(param, 'superpixel', 2500);
compact = getParameter(param, 'compactness', 10);

img8 = imrescale(img, 'class', 'uint8');

[imglab, ~] = slicmex(img8, nspxls, compact);

end



