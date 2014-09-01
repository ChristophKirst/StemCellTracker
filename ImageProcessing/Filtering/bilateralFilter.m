function out = bilateralFilter(image, hsize, sigma_space, sigma_intensity, padding)
%
% out = bilateralFilter(image, hsize sigma_space, sigma_intensity, padding)
%
% description:
%     applies bilateral gaussian filter on image
%
% input:
%     image            image to filter
%     hsize            h x w (x l) filter size
%     sigma_space      spatial kernel width
%     sigma_intensity  intensity kernel width
%     padding          padding method   
%
% output:
%     out              filtered image
%
% reference:
%     C. Tomasi and R. Manduchi. Bilateral Filtering for Gray and Color Images. 
%     In Proceedings of the IEEE International Conference on Computer Vision, 1998.
%
% See also: meanShiftFilter, linearFilter

% initialize

dim = ndims(image);

if dim < 2 || dim > 4
   error('bilateralFilter: expect 2d or 3d gray scale image');
end

image = double(image);

if dim ==2
   [n,m] = size(image);
else
   [n,m,k] = size(image);
end

if length(hsize) < dim
   hsize = repmat(hsize(1), dim, 1);
else
   hsize = hsize(1:dim);
end
[ol, or] = filteroffsets(hsize);
mo = max([ol;or]);

if nargin < 3
   sigma_space = hsize /(2 * sqrt(2 * log(2)));
end
if length(sigma_space) < dim
   s = repmat(sigma_space(1), dim, 1);
else
   s = sigma_space(1:dim);
end
s2 = 2 * s.^2;

if nargin < 4 
   % some heuristic guess 1.1 standard deviation
   si  = 1.1 * std(image(:));  
else
   si = sigma_intensity(1);
end
si2 = 2*si^2;

if nargin < 5
   padding = 'replicate';
end

% filter
if dim == 2

   out = zeros(n,m);       
   blw = zeros(n,m);    

   imgpad = padarray(image, mo, padding);

   for i1 = -ol(1):or(1)  % loop over shifts                        
     for i2 = -ol(2):or(2)       
         ws = exp(-sum([i1 i2].^2./s2));                
         ms  = imgpad( (1+mo(1)+i1):(n+mo(1)+i1), (1+mo(2)+i2):(n+mo(2)+i2) );   % shift image
         wr  = exp(-((image-ms).^2)/si2);      % filter in intensity
         wr  = ws*wr;                        
         out = out + (wr.*ms) ;              % accumulate bilateral summation
         blw = blw + wr ;                    % accumulate bilateral weighting
     end
   end
   out =out./blw;                            
   out = out / max(out(:));

else % dim == 3

   out = zeros(n,m,p);       
   blw = zeros(n,m,p);    

   imgpad = padarray(image, mo, padding);

   for i1 = -ol(1):or(1)  % loop over shifts                        
     for i2 = -ol(2):or(2)
        for i3 = -ol(3):or(3)
           ws = exp(-sum([i1 i2 i3].^2./s2));    
           ms  = imgpad( (1+mo(1)+i1):(n+mo(1)+i1), (1+mo(2)+i2):(m+mo(2)+i2), (1+mo(3)+i3):(p+mo(3)+i3) );   % shift image
           wr  = exp(-((image-ms).^2)/si2);      % filter in intensity
           wr  = ws*wr;                        
           out = out + (wr.*ms) ;              % accumulate bilateral summation
           blw = blw + wr ;                    % accumulate bilateral weighting
       end
     end
   end

   out =out./blw;                            
   out = out / max(out(:));
end