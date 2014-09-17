function troughs = findTroughs(x, delta, k)
%
% troughs = findTroughs(X)
%
% description:
%    finds troughs in time series x
%
% input:
%    x        time series
%    delta    min relative change to detect trough
%    k        returns at most k troughs
%
% output:
%    troughs  positions of troughs
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

troughs = []; s = [];
i = 0; d = 0; xu = x(1); xd = x(1);

while i < n
   i = i + 1; xi = x(i);
      switch d
         case 1
            if xi >= xu % still growing
               xu = xi;
            elseif  xi <= xu - delta  % strong enough decrease
               xd = xi; d = -1; s= i;
            end      
         case -1
            if xi <= xd % still decaying
               s= i; xd = xi;
            elseif xi >= xd + delta % strong enough growth
               troughs = [troughs, s];  %#ok<AGROW>
               xu = xi; 
               d = 1;
               if length(troughs) >= k 
                  return
               end
            elseif xi == xd
               s = [s, i]; %#ok<AGROW>
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



