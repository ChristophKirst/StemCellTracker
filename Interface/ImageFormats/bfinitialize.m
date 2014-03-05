function bfinitialize(bfpath)
%
% [status, version] = bfinitialize(bfpath)
%
% description: 
%    initalizes loci_tool.jar for the use with matlab searching in bfpath
%
% input:
%    bfpath    (optional) a hint to the loci_tools path
%
% output:
%    status    true if loci_tools.jar is in the Java class path.
%    version   version of Bio-Formats if loci_tools.jar, '' otherwise
%
% See also: ijinitialize, ijpath, bfpath

% note: modifed from bfmatlab

if javacheckclasspath('loci_tools.jar')
   status = true;
end

% try to find loci_tools.jar in the specified path or all its subfolders
if ~status && nargin > 0
   loci = findloci(bfpath);
   if ~isempty(loci)
      javaaddjar(loci);
      status = true;
   end
end

% automatic detection: check this path
if ~status
    loci = fullfile(fileparts(mfilename('fullpath')), 'loci_tools.jar');
    if exist(loci, 'file')
      javaaddjar(loci);
      status = true;
    end
end

% still not found, search in imagej folder
if ~status
    bfpath = ijpath;
    if ~isempty(bfpath)
       loci = findloci(bfpath);
       if ~isempty(loci)
         javaaddjar(loci);
         status = true;
       end
    end
end

if ~status
   error('bfinitialize: could not find loci_tools.jar!')
else
   v = eval('loci.formats.FormatTools.VERSION');
   fprintf('bfinitialize: loci_tools.jar installed: version: %s\n', char(v))
end

end


function loci = findloci(bfpath)
   if bfpath(end) ~= filesep
      bfpath = [bfpath filesep];
   end
   bfpath = fileparts(bfpath);

   loci = [bfpath 'loci_tools.jar'];
   
   if exist(lfile, 'file')
      return
   end
   
   d = dir(bfpath);
   for i=1:length(d)
      if d(i).isdir
         loci = findloci(fullfile(bfpath, d(i).name));
      end
      if ~isempty(loci)
         return
      end
   end
   
   loci = '';
end

