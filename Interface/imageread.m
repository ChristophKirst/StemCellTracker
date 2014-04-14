function img = imageRead(cmdformat, imagespec)
%
% img = imageRead(cmdformat, imagespec)
%
% description:
%    opens an image specified by taking the comand format cmdforamt
%    and replacing it with the image specifications imagespec
%    replacements are done as follows:
%    of a field with name xxx in the struct imagespec the substring 
%       <xxx>    is replaced by the value, 
%       <xxx,k>  is replaced by the value xxx using k digits with trailling zeros
%       <#n>     is replaced by the nth element of cell imagespec
%       <#n,k>   is replaced by the value xxx using k digits with trailling zeros
%
% input:
%    cmdformat     string of the command that opens the image
%    imagespec     struct wiht the specifications
%
% See also: imageReadCommand num2str0

cmd = imageReadCommand(cmdformat, imagespec);
img = eval(cmd);

end