function [ids, nids] = imfrmtParseRangeToIndex(si, frmt, varargin)
%
% [ids, nids] = imfrmtParseRangeToIndex(si, frmt, param)
%
% description: 
%    parses parameter in varargin keeps and reformat paramter as in format frmt / size si
%    uses last format specification if multiple
%    works on a single size and fmrt only !
%
% input:
%    si    size
%    frmt  format
%    param parameter
%
% output:
%    ids   indices 
%    nids  number of indices
%
% note:
%    used for the bio formats interace to parse range specs


param = parseParameter(varargin);

ids = getParameter(param, upper(frmt), []);
if isempty(ids)
   ids = getParameter(param, lower(frmt), []);
   if ~isempty(ids)
      ids = si - ids + 1;
   end
end
if isempty(ids)
   ids = 1:si;
end
if iscell(ids)
   ids = cell2mat(ids);
end

if nargout > 1
   nids = length(ids);
end

end