function ilstart()
%
% ilstart(clfile)
%
% description:
%    starts the ilastik software to generate a classifier
%

[ipath, ilrun] = ilpath();

if isunix() 
   system([fullfile(ipath, ilrun) '&'])
else
   system(fullfile(ipath, ilrun))
end