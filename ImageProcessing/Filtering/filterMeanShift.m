function out = filterMeanShift(image, ksize, intensity_width, iterations, padding)
%
% out = filterMeanShift( image, ksize, intensity_width, iterations, padding )
%
% description:
%    applies mean shift filter to image
%
% input:
%    image              image to be filtered
%    ksize              h x w (x l) size of the filter
%    intensity_width    intensity window size
%    iterations         number of repititions
%  
% output:
%     out               filtered image
%
% See also: filterBilateral, filterMedian

% initialize

dim = ndims(image);

if dim < 2 || dim > 4
   error('filterMeanShift: expect 2d or 3d gray scale image!');
end

image = double(image);

if dim == 2
   [n,m] = size(image);
else
   [n,m,p] = size(image);
end


if nargin < 2
   ksize = 3;
end
if length(ksize) < dim
   ksize = repmat(ksize(1), dim,1);
else
   ksize = ksize(1:dim);
end
[ol, or] = filteroffsets(ksize);
mo = max([ol; or]);

if nargin < 3
   intensity_width = std(image(:));
end

if nargin < 4
   iterations = 1;
end
if nargin < 5
   padding = 'replicate';
end


out = image;


if dim ==2 

   for it = 1:iterations
      % padding
      imgpad = padarray(out, mo, padding);

      out = zeros(n,m);       
      blw = zeros(n,m);     

      for i1 = -ol(1):or(1)  % loop over shifts                        
        for i2 = -ol(2):or(2)
          %ws = sum([i1 i2].^2/r2);                       % Euclidean distance    
          %if ws < 1                                      % might want to skip and take pixel distance
            %ws = exp(-ws/(s^2)); 
            ms  = imgpad( (1 + mo(1) + i1):(n + mo(1) + i1), (1 + mo(2) + i2):(m + mo(2) + i2) );   % shift image
            wr  = abs(image-ms) < intensity_width;            % filter in intensity                    
            out = out + (wr.*ms) ;              % accumulate bilateral summation
            blw = blw + wr ;                    % accumulate bilateral weighting
          %end
        end
      end
      out=out./blw;                            
      out = out / max(out(:));
   end

else % dim = 3

   for it = 1:iterations
      % padding
      imgpad = padarray(out, mo, padding);

      out = zeros(n,m,p);       
      blw = zeros(n,m,p);     

      for i1 = -ol(1):or(1)  % loop over shifts                        
        for i2 = -ol(2):or(2)
           for i3 = -ol(3):or(3)
             %ws = sum([i1 i2 i3].^2/r2);                    % Euclidean distance    
             %if ws < 1                                      % might want to skip and take pixel distance
               %ws = exp(-ws/(s^2)); 
               ms  = imgpad( (1+mo(1)+i1):(n+mo(1)+i1), (1+mo(2)+i2):(m+mo(2)+i2), (1+mo(3)+i3):(p+mo(3)+i3) );   % shift image
               wr  = abs(image-ms) < intensity_width;            % filter in intensity                    
               out = out + (wr.*ms) ;              % accumulate bilateral summation
               blw = blw + wr ;                    % accumulate bilateral weighting
             %end
           end
        end
      end
      out=out./blw;                            
      out = out / max(out(:));
   end

end


