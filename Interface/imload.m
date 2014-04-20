function img = imload(cmd, fname, param)
%
% img = imload(fname)
% img = imload(fname, param)
% img = imload(cmd, fname, param)
%
% description: 
%    returns an image using the file format fname, command cmd and parameter param
%
% input:
%     cmd     tagged command to read the image, if not specified 'imread_bf(<file>)' is used (only tag is <file>)
%     fname   tagged filename format e.g. 'img_Z<zslice,2>_T<time,2>.tif' (coordinate tags or arbitrary tags in which case the tagnames speficied in param are used)
%     param   (optional)  parameter struct with coordinate entries of the form 
%                         {min, max}, {min, max, delta} or absolute indices as array [i1, i2, ...}
%                         coordinates are
%             .p          p = first spatial coordinate
%             .q          q = second spatial coordinate 
%             .l          l = third spatial coordinate
%             .t          t = time coordinate
%             .c          c = channel coordinate
%             .tagnames   a sturct that maps coordinates names to tagnames in the file
%
% output:
%     img     the image as numeric array \


%% bla






  