function ker = fspecial2(type, ksize, varargin)
%
% ker = fspecial2(type, ksize, varargin)
%
% description:
%      creates sepcial 2d filter kernel in consistent way
%
% input:
%      type      'average', 'gaussian', 'disk', 'laplacian', 'log', or 'dog'
%      ksize     h x w fitler box 
%      varargin  std for gaussian and log filter
% 
% output:
%      ker       2D fitler kernel matrix

if ~ischar(type)
   error('fspecial2: type argument must be an array of characters!')
end
type = lower(type);

if nargin < 2  || isempty(ksize)
   ksize = 3;
end
if length(ksize) < 2
   ksize = [ksize(1) ksize(1)];
else
   ksize = ksize(1:2);
end

[ol,or] = filteroffsets(ksize);
mo = max([ol; or]);

switch type
    
   case 'average'
      
      ker = ones(ksize)/prod(ksize);
        
   case 'gaussian'         
       
      if nargin < 3  || isempty(varargin{1})
         sigma = ksize / 2 / sqrt(2 * log(2));
      else
         sigma = varargin{1};
      end
      if length(sigma) < 2
         sigma = [sigma(1) sigma(1)];
      else
         sigma = sigma(1:2);
      end

      [x,y] = ndgrid(-ol(1):or(1),-ol(2):or(2));
      add = mod(ksize + 1,2)/2;
      x = x + add(1); y= y + add(2);
      
      ker = exp(-(x.*x/2/sigma(1)^2 + y.*y/2/sigma(2)^2));
      ker = ker/sum(ker(:));
  
   case 'sphere'

      add = mod(ksize + 1,2)/2;
      if nargin < 3  || isempty(varargin{1})
         radius = mo - add + 0.5;
      else
         radius = varargin{1};
      end
      if length(radius) < 2
         radius = [radius(1) radius(1)];
      else
         radius = radius(1:2);
      end
      
      [x,y] = ndgrid(-ol(1):or(1),-ol(2):or(2));
      x = x + add(1); y= y + add(2);
      
      ker = 1 - (x.^2 / radius(1)^2 + y.^2 / radius(2)^2);
      ker(ker < 0) = 0;
      ker = ker / sum(ker(:));

   case 'disk'
       
      if nargin < 3  || isempty(varargin{1})
         ring_width = 0;
      else
         ring_width = varargin{1};
      end
      if length(ring_width) < 2
         ring_width = [ring_width(1) ring_width(1)];
      else
         ring_width = ring_width(1:2);
      end

      if nargin < 4  || isempty(varargin{2})
         w_disk = 1;
      else
         w_disk = varargin{2};
      end
      if nargin < 5  || isempty(varargin{3})
         w_ring = -1;
      else
         w_ring = varargin{3};
      end

      add = mod(ksize + 1,2)/2;
      radius_out = max(1.5, mo + add) + 10 * eps;
      radius_in = max(1.5, mo + add - ring_width) + 10 * eps;
      
      [x,y] = ndgrid(-ol(1):or(1),-ol(2):or(2));
      x = x + add(1); y = y + add(2);      
      
      ker = 1 - (x.^2 / radius_out(1)^2 + y.^2 / radius_out(2)^2);
      ker(ker < 0) = 0;
      ker(ker > 0) = w_ring;

      ker2 = 1 - (x.^2 / radius_in(1)^2 + y.^2 / radius_in(2)^2);
      ker2(ker2 < 0) = 0;
      ker2(ker2 > 0) = w_disk - w_ring;

      ker = ker + ker2;


      if isequal(ring_width, [0, 0])
         ker = ker / sum(ker(:));
      end
    
   case 'laplacian'
      
      if nargin < 3  || isempty(varargin{1})
         alpha = 0;
      else
         alpha = varargin{1};
      end
      
      ker = fspecial('laplacian',alpha);

      
   case 'log'
       
       if nargin < 3 || isempty(varargin{1})
          sigma = ksize / 4 / sqrt(2 * log(2));
       else
          sigma = varargin{1};
       end
       if length(sigma) < 2
          sigma = [sigma(1) sigma(1)];
       else
          sigma = sigma(1:2);
       end

       [x,y] = ndgrid(-ol(1):or(1),-ol(2):or(2));
       add = mod(ksize + 1,2)/2;
       x = x + add(1); y= y + add(2);
       
       
       ker = exp(-(x.*x/2/sigma(1)^2 + y.*y/2/sigma(2)^2));
       ker = ker/sum(ker(:));
       arg = (x.*x/sigma(1)^4 + y.*y/sigma(2)^4 -(1/sigma(1)^2 + 1/sigma(2)^2));
       ker = arg.*ker;
       ker = ker - sum(ker(:))/numel(ker);
             
   case 'dog'

      if nargin < 4 || isempty(varargin{2})
         sigma_out = ksize / 2 / sqrt(2 * log(2));
      else
         sigma_out = varargin{2};
      end
      if length(sigma_out) < 2
         sigma_out = [sigma_out(1) sigma_out(1)];
      else
         sigma_out = sigma_out(1:2);
      end

      if nargin < 3 || isempty( varargin{1})
         sigma_in = sigma_out / 1.5;
      else
         sigma_in = varargin{1};
      end
      if length(sigma_in) < 2
         sigma_in = [sigma_in(1) sigma_in(1)];
      else
         sigma_in = sigma_in(1:2);
      end
      
      [x,y] = ndgrid(-ol(1):or(1),-ol(2):or(2));
      add = mod(ksize + 1,2)/2;
      x = x + add(1); y= y + add(2);
      
      ker = exp(-(x.*x/2/sigma_in(1)^2 + y.*y/2/sigma_in(2)^2));
      ker = ker/sum(ker(:));
      
      sub = exp(-(x.*x/2/sigma_out(1)^2 + y.*y/2/sigma_out(2)^2));
      ker = ker - sub / sum(sub(:));
      
      
      
   case 'cone'  % cone
     
      if nargin < 3  || isempty(varargin{1})
         wout = 1;
      else
         wout = varargin{1};
      end
      
      if nargin < 4 || isempty(varargin{2})
         win = 0;
      else
         win = varargin{2};
      end
      
      
      add = mod(ksize + 1,2)/2;
      radius = max(1.5, mo + add) + 10 * eps;
      
      [x,y] = ndgrid(-ol(1):or(1),-ol(2):or(2));
      x = x + add(1); y = y + add(2);      
         
      ker = (1 - sqrt((x.^2 / radius(1)^2 + y.^2 / radius(2)^2)));
      idx = ker <= 0;
      ker = ker * (wout - win) + win;
      ker(idx) = 0;
      
      if win >= 0 && wout >= 0
         ker = ker/sum(ker(:));
      end
      
      
   otherwise
        error('fspecial2: unknown filter type.')
        
end
