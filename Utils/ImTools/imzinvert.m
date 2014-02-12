function stack = imzinvert(stackin)
%
% stack = imzinvert(stack)
%
% description:
%    inverst the z axis in the stack data
%
% input:
%    stackin    input stack
%
% output:
%    stack      stack with inverted z direction
%

stack = stackin;
stack = flipdim(stackin, imdim(stack, 'z'));

end


