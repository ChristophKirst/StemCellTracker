%% Test 3D segmentation with new data



stack = imreadstack('./Test/Images/Stack2/*c001.tif', 'PixelRegion', {[x0 x0+wh], [y0 y0+wh]});

ijplot3d(stack)