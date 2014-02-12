function peaks = findPeaks(x, delta, k)
%
% peaks = findPeaks(X)
%
% description:
%    finds peaks in yime series x
%
% input:
%    x      time series
%    delta  min relative change to detect peak
%    k      returns at most k peaks
%
% output:
%    peaks  positions of peaks
%
% reference: 
%    Todd, Andrews

if nargin < 2
   delta = 1;
end
if nargin < 3
   k = Inf;
end

n = length(x);

peaks = []; s = [];
i = 0; d = 0; xu = x(1); xd = x(1);

while i < n
   i = i + 1; xi = x(i);
      switch d
         case 1
            if xi > xu % still growing
               s = i; xu = xi;
            elseif  xi <= xu - delta  % strong enough decrease
               peaks = [peaks, s]; xd = xi; d = -1; %#ok<AGROW>
               if length(peaks) >= k 
                  return
               end
            elseif xi == xu % plateau 
               s = [s, i]; %#ok<AGROW>
            end
         case -1
            if xi <= xd % still decaying
               xd = xi;
            elseif xi >= xd + delta % strong enough growth
               xu = xi; d = 1; s= i;
            end
         case 0 % for initial state
            if xi <= xu - delta % trend is decaying
               d = -1; s = i;
            elseif xi >= xd + delta; % trend is rising
               d = 1;
            end
            if xi < xd  % update max / min
               xd = xi;
            elseif xi > xu
               xu = xi;
            end
      end
end

end


% clipboard from/ converted partially from Mathematica code
% 
% function [peaks troughs] = findPeakTroughs(x, delta)
% 
%    n = length(x);
%    p = []; s = []; t = {};
%    i = 0; d = 0; qa = x(1); qb = x(1);
%    while i < n,
%       i = i+1; qi = x(i);
%       switch d
%          case 1
%             if qa < qi
%                s = {i}; qa = qi;
%             elseif qa == qi
%                s = {s, i};
%             elseif qa >= qi + delta
%                p = [p, s]; s = i; qb = qi; d = -1
%             end
%          case -1
%             if qi <=  qb
%                s = i; qb = qi;
%             elseif qi == qb
%                s = [s, i];
%             elseif qi >= qb + delta
%                t = [t, s]; s = i; qa = qi; d = 1;
%             end
%          case 0
%             if qa >= qi + delta
%                d = -1;
%             elseif qi >= qb + delta
%                d = 1;
%             end;
%             if qa < qi
%                qa = qi;
%             elseif qi < qb
%                qb = qi;
%             end
%             s = i;
%       end
%    end
%    
%    peaks = p;
%    troughs = t;
% end
% 
% 
% 


% 
% (* detects the amplyidues between peak and succeding trough *)
% PeakAmplitudes[db_, pks0_, trs0_] := Module[{amp, pks = pks0, trs = trs0},
%    If[Length[pks] == 0, Return[{}]];
%    If[Length[trs] > 0 && trs[[1]] < pks[[1]], trs =Drop[trs,1]];
%    If[Length[trs] < Length[pks],
%       trs = Append[trs, pks[[-1]] - 1 + Position[db[[pks[[-1]] ;; -1]],  Min[db[[pks[[-1]] ;; -1]]], 1, 1][[1, 1]]]
%    ];
%    Return[db[[pks]] - db[[trs]]]
% ];
% 
% 
% TroughAmplitudes[db_, pks0_, trs0_] := Module[{amp, pks = pks0, trs = trs0},
%    If[Length[trs] == 0, Return[{}]];
%    If[Length[pks] > 0 && pks[[1]] < trs[[1]], pks =Drop[pks,1]];
%    If[Length[pks] < Length[trs],
%       pks = Append[pks, trs[[-1]] - 1 + Position[db[[trs[[-1]] ;; -1]],  Max[db[[trs[[-1]] ;; -1]]], 1, 1][[1, 1]]]];
%    Return[db[[pks]] - db[[trs]]]
% ];
