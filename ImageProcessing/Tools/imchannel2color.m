function c = imchannel2color(name)


switch lower(name)
   case 'dapi' 
      c = 'b';
   case 'gfp'
      c = 'g';
   case 'rfo'
      c = 'r';
   case 'yfp'
      c = 'y';

   otherwise
      c = 'gray';
end

end