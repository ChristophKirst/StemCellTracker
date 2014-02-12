function mask = findObjectCenter(image, param)
%
% mask = findObjectCenter(image, param)
%
% description:
%    locates possible centers of objects in an image
%
% input:
%    image   the image to locate object centers
%    param   parameter struct with
%            .method         method for finding the center
%            .parameter      parameter structure for sub-method
%
% output:
%    mask    a bw image containing possible centers as white pixels
%

method = getParameter(param, {'method'}, 'log');


switch method
   case 'log'
      filter = logFilter(image, param.parameter);
      
   case '
      
      
   otherwise
      error(['findObjectCenter: unknown method: ' method '!']);
end

end



