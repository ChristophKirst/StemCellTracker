function stack = imzreverse(stackin)
%
% stack = imzreverse(stack)
%
% description:
%    reverses the z axis in the stack data
%
% input:
%    stackin    input stack
%
% output:
%    stack      stack with inverted z direction
%

stack = flip(stackin, imdim(stackin, 'z'));

end


