function ireader = imreadBFReader(name, varargin)
%
% ireader = imreadBFReader(name, varargin)
%
% description:
%    parses input and returns a java bio-formats image reader class
%
% input: 
%    name    (optional) filename
%
% output:
%    ireader image reader

% Prompt for a file if no input
if nargin == 0
  [file, path] = uigetfile(bffileextensions, 'Choose a file to open');
  name = [path file];
end

if ischar(name)  
   
   % load via file name
   %if ~exist(name, 'file') && ~isdir(name)
   %   error('imreadBFReader: file does not exists: %s', name)
   %end
   
   try 
      % initialize logging
      % loci.common.DebugTools.enableLogging('INFO');

      % Get the channel filler
      ireader = loci.formats.ChannelSeparator(loci.formats.ChannelFiller());
   
      % initialize file stitcher
      % ireader = loci.formats.FileStitcher(ireader);
      % ireader.setGroupFiles(false);
   
      % initialize meta data
      OMEXMLService = loci.formats.services.OMEXMLServiceImpl();
      ireader.setMetadataStore(OMEXMLService.createOMEXMLMetadata());
      %ireader.setMetadataStore(loci.formats.MetadataTools.createOMEXMLMetadata());
   
      % initialize reader from file name
      ireader.setId(name);
   catch
      error('imreadBFReader: cannot open image file %s via bio-formats, try to run bfinitialize first!', name);
   end
  
elseif isa(name, 'loci.formats.ChannelSeparator') 
   
   ireader = name;
   
end

end
