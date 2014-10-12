function img = imreadAVI(fname, varargin)
%
% img = imreadAVI(fname, varargin)
%
% description:
%     reads an AVI file as XYT array

param = parseParameter(varargin{:});




trange = getParameter(param, 'time', []);

vid = VideoReader(fname);
vidHeight = vid.Height;
vidWidth = vid.Width;
vidn = vid.NumberOfFrames;
vidn = 1:vidn;
if ~isempty(trange)
   vidn(vidn < trange(1)) = [];
   vidn(vidn > trange(2)) = [];
end

img = zeros(vidHeight, vidWidth, length(vidn), 3);

k = 1;
for t = vidn
    img(:,:,k,:)  = read(vid,t);
    k = k + 1;
end


end

