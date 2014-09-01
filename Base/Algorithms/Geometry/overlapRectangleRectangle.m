function a = overlapRectangleRectangle(a,b)
   a1 = a(1,:); a2 = a(2,:);
   b1 = b(1,:); b2 = b(2,:);
   
   ov = zeros(2,length(a1));
   ov(1,:) = max(a1,b1);
   ov(2,:) = min(a2,b2);
   
   if any(ov(2,:)-ov(1,:) < 0)
      a = 0;
   else
      a = prod(ov);
   end
end