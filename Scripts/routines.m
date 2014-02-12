%% assembly of crucial image procesing steps




%% smoothing 
                sigma = 1;
                FiltLength = 2*sigma;
                [x,y] = meshgrid(-FiltLength:FiltLength,-FiltLength:FiltLength);   % Filter kernel grid
                f = exp(-(x.^2+y.^2)/(2*sigma^2));f = f/sum(f(:));                 % Gaussian filter kernel
                %                BlurredImage = conv2(OrigImage,f,'same');                             % Blur original image
                %%% This adjustment prevents the outer borders of the image from being
                %%% darker (due to padding with zeros), which causes some objects on the
                %%% edge of the image to not  be identified all the way to the edge of the
                %%% image and therefore not be thrown out properly.
                BlurredImage = conv2(OrigImage,f,'same') ./ conv2(ones(size(OrigImage)),f,'same');
                
                






