function ker = fspecial3(type, ksize, varargin)
%
% ker = fspecial3(type, ksize, varargin)
%
% description:
%      creates sepcial 3d filter kernel in consistent way
%
% usage:
%       fspecial3('average', ksize)
%       fspecial3('disk', ksize, ring_width, w_disk, w_ring)
%
% input:
%      type      'average', 'gaussian', 'sphere', 'disk', 'laplacian', 'log', or 'dog'
%      ksize     h x w x l fitler box 
%      varargin  std for gaussian and log filter
% 
% output:
%      ker       3D fitler kernel matrix

if ~ischar(type)
   error('fspecial3: type argument must be an array of characters!')
end
type = lower(type);

if nargin < 2
   ksize = 3;
end
if length(ksize) < 3
   ksize = [ksize(1) ksize(1), ksize(1)];
else
   ksize = ksize(1:3);
end

[ol,or] = filteroffsets(ksize);
mo = max([ol; or]);

switch type
    
   case 'average'
      
      ker = ones(ksize)/prod(ksize);
        
   case 'gaussian'         
       
      if nargin < 3
         sigma = ksize / 2 / sqrt(2 * log(2));
      else
         sigma = varargin{1};
      end
      if length(sigma) < 3
         sigma = [sigma(1) sigma(1) sigma(1)];
      else
         sigma = sigma(1:3);
      end

      [x,y,z] = ndgrid(-ol(1):or(1),-ol(2):or(2), -ol(3):or(3));
      ker = exp(-(x.*x/2/sigma(1)^2 + y.*y/2/sigma(2)^2) + z.*z/2/sigma(3)^2);
      ker = ker/sum(ker(:));
  
   case 'sphere'
      
      if nargin < 3
         radius = mo;
      else
         radius = varargin{1};
      end
      
      if length(radius) < 3
         radius = [radius(1) radius(1) radius(1)];
      else
         radius = radius(1:3);
      end
      
      [x,y,z] = ndgrid(-ol(1):or(1),-ol(2):or(2),-ol(3):or(3));
      x = x+0.5; y = y + 0.5; z = z + 0.5;

      ker = 1 - (x.^2 / radius(1)^2 + y.^2 / radius(2)^2 + z.^2 / radius(3)^2);
      ker(ker < 0) = 0;
      ker = ker / sum(ker(:));

   case 'disk'
       
      if nargin < 3
         ring_width = 0;
      else
         ring_width = varargin{1};
      end
      if length(ring_width) < 3
         ring_width = [ring_width(1) ring_width(1) ring_width(1)];
      else
         ring_width = ring_width(1:3);
      end

      if nargin < 4
         w_disk = 1;
      else
         w_disk = varargin{2};
      end
      if nargin < 5
         w_ring = -1;
      else
         w_ring = varargin{3};
      end

      radius_in = max(0, mo - ring_width);
      radius_out = min(mo, radius_in + ring_width);

      [x,y,z] = ndgrid(-ol(1):or(1),-ol(2):or(2),-ol(3):or(3));
      x = x+0.5; y = y + 0.5; z = z + 0.5;

      ker = 1 - (x.^2 / radius_out(1)^2 + y.^2 / radius_out(2)^2 + z.^2 / radius_out(3)^2);
      ker(ker < 0) = 0;
      ker(ker > 0) = w_ring;

      ker2 = 1 - (x.^2 / radius_in(1)^2 + y.^2 / radius_in(2)^2 + z.^2 / radius_in(3)^2);
      ker2(ker2 < 0) = 0;
      ker2(ker2 > 0) = w_disk - w_ring;

      ker = ker + ker2;

      if (w_ring == 0) || isequal(ring_width, [0, 0, 0])
         ker = ker / sum(ker(:));
      end
    
   case 'laplacian'
      
      %if nargin < 3
      %   alpha = 0;
      %else
      %   alpha = varargin{1};
      %end
      
      ker = zeros(3,3,3);
      ker(:,:,1) = [0 3 0;3 10 3;0 3 0];
      ker(:,:,3) = ker(:,:,1);
      ker(:,:,2) = [3 10 3;10 -96 10;3 10 3];
     
   case 'log'
       
      if nargin < 3
         sigma = ksize / 4 / sqrt(2 * log(2));
      else
         sigma = varargin{1};
      end
      if length(sigma) < 3
         sigma = [sigma(1) sigma(1) sigma(1)];
      else
         sigma = sigma(1:3);
      end

      [x,y,z] = ndgrid(-ol(1):or(1),-ol(2):or(2), -ol(3):or(3));
      ker = exp(-(x.*x/2/sigma(1)^2 + y.*y/2/sigma(2)^2 + z.*z/2/sigma(3)^2));
      ker = ker/sum(ker(:));
      arg = (x.*x/sigma(1)^4 + y.*y/sigma(2)^4 + z.*z/sigma(3)^4 - (1/sigma(1)^2 + 1/sigma(2)^2 +1/sigma(3)^2));
      ker = arg.*ker;
      ker = ker - sum(ker(:))/numel(ker);
             
   case 'dog'

      if nargin < 4 || isempty(varargin{2})
         sigma_out = ksize / 2 / sqrt(2 * log(2));
      else
         sigma_out = varargin{2};
      end
      if length(sigma_out) < 3
         sigma_out = [sigma_out(1) sigma_out(1) sigma_out(1)];
      else
         sigma_out = sigma_out(1:3);
      end

      if nargin < 3 || isempty( varargin{1})
         sigma_in = sigma_out / 1.5;
      else
         sigma_in = varargin{1};
      end
      if length(sigma_in) < 3
         sigma_in = [sigma_in(1) sigma_in(1) sigma_in(1)];
      else
         sigma_in = sigma_in(1:3);
      end
      
      [x,y,z] = ndgrid(-ol(1):or(1),-ol(2):or(2),-ol(3):or(3));
      ker = exp(-(x.*x/2/sigma_in(1)^2 + y.*y/2/sigma_in(2)^2 + z.*z/2/sigma_in(3)^2));
      ker = ker/sum(ker(:));
      
      sub = exp(-(x.*x/2/sigma_out(1)^2 + y.*y/2/sigma_out(2)^2 + z.*z/2/sigma_out(3)^2));
      ker = ker - sub / sum(sub(:));
      
   otherwise
        error('fspecial3: unknown filter type.')
        
end
