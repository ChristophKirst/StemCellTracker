function plotMatchedObjects(match)
%
% plotMatchedObjects(match)
%
% description:
%    plot result of matching
%
% input:
%    match     Match object
% 
% See also: Match

if ~isa(match, 'Match')
   error('expect match to be of class Match')
end


xyz0 = [ match.objects0.r ];
s0 = [ match.objects0.volume ];
s0 = s0 / max(s0) .* 200;

xyz1 =  [ match.objects1.r ];
s1 = [ match.objects1.volume ];
s1 = s1 / max(s1) .* 200;

hold on
grid on

if match.dim == 2
   h0 = scatter(xyz0(1,:),xyz0(2,:), s0, [1 0 0]);
   h1 = scatter(xyz1(1,:),xyz1(2,:), s1, [0 0 1]);
   
   for i=1:match.n0
      j = match.match(i);
      if j>0
         hl = line([xyz0(1,i) xyz1(1,j)], [xyz0(2,i) xyz1(2,j)]);
      end
   end
else
   h0 = scatter3(xyz0(1,:),xyz0(2,:),xyz0(3,:), s0, [1 0 0]);
   h1 = scatter3(xyz1(1,:),xyz1(2,:),xyz1(3,:), s1, [0 0 1]);
   
   for i=1:match.n0
      j = match.match(i);
      if j>0
         hl = line([xyz0(1,i) xyz1(1,j)], [xyz0(2,i) xyz1(2,j)], [xyz0(3,i) xyz1(3,j)]);
      end
   end
end


% some legends 
xlabel('x'); ylabel('y'); zlabel('z');
legend([h0, h1, hl], 'r0', 'r1', 'match');

hold off


end

