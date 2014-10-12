function stack = imread_hdr(info)
%
% stack = imread_hdr(info)
%
% description:
%     function for reading volume of HDR/IMG Analyze ( .hdr ) volume file
% 
% input: 
%     info    filename or file info
%
% See also: analyze75info analyze75read

if (~isstruct(info))
   info=analyze75info(info); 
end

stack = analyze75read(info);
stack = flip(stack, 1);

end
