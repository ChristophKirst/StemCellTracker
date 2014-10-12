function moviename = movieFromImages(imgs, moviename, varargin)
%
% moviename = movieFromImages(imgs, moviename)
%
% description: 
%     reads images from the list and writes them to a movie
%
% input:
%     imgs         list of images, filenamesor an ImageSource
%     moviename    file name of the movie
%     param        parameter struct with entries
%                  .type       movie type 'avi' or 'mp4' or as in VideoWriter
%                  .framerate  the frame rate in frames/second
%                  .norm       divide by this to get image in range 0 - 1 
%
% output:
%     moviename    file name of the movie
%
% See also: VideoWriter

param = parseParameter(varargin{:});

type = getParameter(param, 'type', 'avi');
switch type
   case 'avi' 
      type = 'Motion JPEG AVI';
%    case 'mp4'
%       %type = 'MPEG-4';
%       type = 'Motion JPEG 2000';
%       if isempty(param)
%          clear param
%          param.MJ2BitDepth = 8;
%       elseif ~isfield(param, 'MJ2BitDepth')
%          param.MJ2BitDepth = 8;
%       end     
end

fr = getParameter(param, 'framerate', 5);
nrm = getParameter(param, 'norm', []);

vw = VideoWriter(moviename, type);
vw.FrameRate = fr;

if ~isempty(param)
   fn = fieldnames(param);
   for i = 1:length(fn)
      if isprop(vw, fn{i})
         vw.(fn{i}) = param.(fn{i});
      end
   end
end

open(vw);

if iscellstr(imgs) || iscell(imgs)
   nframes = length(imgs);
elseif isa(imgs, 'ImageSource')
   nframes = imgs.cellNumber;
else
   error('movieFromImages: expects cell array or ImageSource as inputs!')
end
   
      
   
  
fprintf(1,'writing image:               ');
for i = 1:nframes
   fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b% 5d / % 5d',i, nframes);
   
   if isa(imgs, 'ImageSource')
      img = imgs.cellData(i);
   elseif iscellstr(imgs);
      img = imreadBF(imgs{i});
   else
      img = imgs{i};
   end
   
   switch class(img)
      case 'double'
         if ~isempty(nrm)
            img = img / nrm;
         end
      otherwise
         if ~isempty(nrm)
            nrm = double(intmax(class(img)));
         end
         img = double(img);
         img = img / nrm;
   end
   img = imclip(img, 0,1);
   writeVideo(vw,img);
   
end
fprintf('\n');
   
end

