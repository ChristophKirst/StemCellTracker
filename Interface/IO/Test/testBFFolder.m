%

clear all
close all
clc
initialize
bfinitialize

%%

% test direcory reader
clc

name = java.util.regex.Pattern.compile('./Test/Images/hESCells_Stack/W1F127T0001Z14C1.tif')

name= ('./Test/Images/hESCells_Stack/W1F127T0001Z01C1.tif')

      % initialize logging
      % loci.common.DebugTools.enableLogging('INFO');

      % Get the channel filler
      % ireader = loci.formats.ChannelSeparator(loci.formats.ChannelFiller());
   
      % initialize file stitcher
      %ireader = loci.formats.FileStitcher(ireader);
      % ireader.setGroupFiles(false);
      
      

      
      %ireader.setId(name);
      
      
      
      %%
      
      stitch = true;
      name= ('/home/ckirst/Science/Simulation/Matlab/StemCell3D/Test/Images/hESCells_<Stack,Test>/W1F127T0001Z<00-14>C<1-2>.tif')
      
      name= ('/home/ckirst/Science/Simulation/Matlab/StemCell3D/Test/Images/hESCells_Folder/f<1-2>_t<01-02>/f<1-2>_t<01-02>_z<01-02>')
      
      name= ('/home/ckirst/Science/Simulation/Matlab/StemCell3D/Test/Images/hESCells_Folder/f1_t1/f1_t01_z01.tif')
      
      name= ('/home/ckirst/Science/Simulation/Matlab/StemCell3D/Test/Images/mESCells_Wnt/wnt_t01_z01.tif')
      id  = name
      
      
      exist(name, 'file')
      
      %name= ('/home/ckirst/Science/Simulation/Matlab/StemCell3D/Test/Images/hESCells_Stack/')
            
      reader = loci.formats.ImageReader()
      
      if (stitch) 
         reader = loci.formats.FileStitcher(reader, true);
         
          f = loci.common.Location(name)

          f.exists()
          
         if (~f.exists()) 

            reader.setUsingPatternIds(true);
            pat = name;

         else 
            
            pat = loci.formats.FilePattern.findPattern(f)
            
         end
         
         if (pat ~= []) 
               id = pat;
         end
      end
   
      reader.setGroupFiles(true)
      
      id
      
      %%
      
      
      reader.setId(id);
      
      
      %%
      fp = reader.getFilePattern
      
      
     
      fp.getPattern

      %%
      
      
      ireader.getFilePattern
      
      %%
      
      % = ireader.findPatterns('./Test/Images/hESCells_Stack/*.tif')
      
      
      
      %      
      
   
      % initialize meta data
      OMEXMLService = loci.formats.services.OMEXMLServiceImpl();
      ireader.setMetadataStore(OMEXMLService.createOMEXMLMetadata());
      %ireader.setMetadataStore(loci.formats.MetadataTools.createOMEXMLMetadata());
   
      % initialize reader from file name
 