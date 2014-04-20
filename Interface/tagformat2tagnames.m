function [tnames, tw] = tagformat2tagnames(tfrmt)
%
% [tnames, tw] = tagformat2tagnames(tfrmt)
%
% description:
%    finds all tag names and possible width specifications in a tagged format string
%
% input:
%    tfrmt   tagged format
% 
% outut:
%    tnames  tagged names
%    tw      width of tag, 0 if not specified.

tgs = regexp(tfrmt, '<(?<name>\w*)\s*(?<k>(,\s*\d*?|\s*?))\s*>', 'names');

tnames = {tgs.name};
tw = zeros(1, length(tnames));

for i = length(tnames):-1:1
   if ~isempty(tgs(i).k)
      tk = tgs(i).k;
      tk(tk == ',') = [];
      tw(i) = str2double(tk);
   end
end

end