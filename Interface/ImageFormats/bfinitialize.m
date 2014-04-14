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

status = false;
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

%still not found, search in bfinitialize's folder
if ~status
   loci = findloci(fileparts(mfilename('fullpath')));
   if ~isempty(loci)
      javaaddjar(loci);
      status = true;
   end
end

%  automatic detection: check imagej's jars path
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

   try 
      v = bfversion();
      fprintf('bfinitialize: loci_tools.jar installed: version: %s\n', v)
   catch
      error('bfinitialize: error while trying to run bf loci_tools.jar')
   end

end

% som classes depend on imagej ! -> matlab gives cannot foind loc class errro with is a non-sense error message
ijinitialize;


end


function loci = findloci(bfpath)
   %if bfpath(end) ~= filesep
   %   bfpath = [bfpath filesep];
   %end
   %bfpath = fileparts(bfpath);

   loci = fullfile(bfpath, 'loci_tools.jar');
   
   if exist(loci, 'file')
      return
   end
   
   d = dir(bfpath);
   
   for i=1:length(d)
      if d(i).isdir && ~any(strcmp(d(i).name, {'.', '..'}))
         loci = findloci(fullfile(bfpath, d(i).name));        
         if ~isempty(loci)
            return
         end
      end

   end
   
   loci = '';
end

