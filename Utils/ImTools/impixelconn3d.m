function ker = impixelconn3d(conn)
%
% ker = impixelconn3d(conn)
%
% description: 
%    returns pixel connectivity
%
% input:
%    conn  the neighbourhood connectivity in 3d: 6, 18, 26
% 
% output:
%    ker   the kernel, ones indicating the connected neighbours

ker = zeros(3,3,3);

switch conn
   case 6
      ker(:,:,1) = [ 0, 0, 0;
                     0, 1, 0;
                     0, 0, 0];
                   
      ker(:,:,2) = [ 0, 1, 0;
                     1, 0, 1;
                     0, 1, 0];
                   
      ker(:,:,3)= ker(:,:,1);
      
   case 18
      ker(:,:,1) = [ 0, 1, 0;
                     1, 1, 1;
                     0, 1, 0];
                   
      ker(:,:,2) = [ 1, 1, 1;
                     1, 0, 1;
                     1, 1, 1];
                   
      ker(:,:,3)= ker(:,:,1); 
      
   case 26
      ker(:,:,1) = [ 1, 1, 1;
                     1, 1, 1;
                     1, 1, 1];
                   
      ker(:,:,2) = [ 1, 1, 1;
                     1, 0, 1;
                     1, 1, 1];
                   
      ker(:,:,3)= ker(:,:,1);
      
   otherwise
      error('impxielconn3d: input not a pixel connectivity');
end
          
      
      
      