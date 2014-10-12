function setPath()
%
% path = setPath()
%
% description:
%     sets all necessary paths
%

includepaths = {
     '/Base/Classes', ...
     '/Base/ImageFormat',...
     '/Base/Parameter',... 
     '/Base/Tags',...
     '/Base/Algorithms',...
     '/Base/Algorithms/Clustering',...    
     '/Base/Algorithms/Geometry',...     
     '/Base/Algorithms/GraphTheory',...     
     '/Base/Algorithms/SignalProcessing',...  
     '/Base/Utils',...
     ...
     '/ImageProcessing/', ...
     '/ImageProcessing/Filtering', ...
     '/ImageProcessing/Photometry', ...
     '/ImageProcessing/Segmentation', ...
     '/ImageProcessing/ShapeDetection', ...
     '/ImageProcessing/Stitching', ...
     '/ImageProcessing/Tracking',...
     '/ImageProcessing/Thresholding', ...
     '/ImageProcessing/Tools',...
     ...
     '/Interface',...
     '/Interface/ImageJ',...
     '/Interface/Imaris',...
     '/Interface/IO',... 
     '/Interface/MetaMorph',...
     '/Interface/Movies',...
     '/Interface/Python',...
     '/Interface/Ilastik',...
     '/Interface/Hugin',...
     ...
     '/Test'...
};
 
basepath = fileparts(mfilename('fullpath'));

for p = 1:length(includepaths)
   includepaths{p} = fullfile(basepath, includepaths{p});
end

addpath(includepaths{:});
     
          
% compability to matlab previous versions
v = version('-release');
if length(v) >= 4
   v = v(1:4);
   if strcmp(v, '2012')
      addpath(fullfile(basepath, '/Base/Utils/External/Matlab2012'));
   end
end

      
end

