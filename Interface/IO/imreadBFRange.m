function rgs = imreadBFRange(name, varargin)
%
% pos = imreadBFRange(name, varargin)
%
% decription:
%    returns internal coordiante ranges for bio-format image
%    checks validity of tag ranges
%    uses last valid tag id in varargin
%
% input:
%    name    filename or reader
%    param   parameter struct with entries
%           .S / .s       series ids 
%           .T / .t       ids of time frames to import [tid1, tid2, ...] or {tid1,tid1, ...} ([] = all)
%           .C / .c       ids of channels to import    [cid1, cid2, ...] or {cid1, cid2,...} ([] = all)
%           .Z / .y       pixel ids in z / l direction [lid1, lid2, ...] or {lid1,lid1, ...} ([] = all)
%           .X / .x       pixel ids in x / p direction [pid1, pid2, ...] or {pid1,pid1, ...} ([] = all)
%           .Y / .y       pixel ids in y / q direction [qid1, qid2, ...] or {qid1,qid1, ...} ([] = all)
%           .format       specifies the format of the ranges to be returned
%             
% output:
%    rgs    cordinate ranges
% 
% note:
%    full X and Y ranges are not returned unless specifed in format !
%

param = parseParameter(varargin);

ireader = imreadBFReader(name);

frmt = getParameter(param, 'format', 'ZCTS');
si = imreadBFSize(ireader, frmt);

n = length(frmt);
fnames= lower(fieldnames(param));

for i = 1:n
   if ~any(ismember(fnames, lower(frmt(i))))
      param.(frmt(i)) = 1:si(i);
   end
end

rgs  = imfrmtParseRange(si, frmt, param);

end
