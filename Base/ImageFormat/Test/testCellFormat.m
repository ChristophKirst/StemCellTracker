
% example: 

cfrmt = 'UV'; 
dfrmt = 'XY';



for i = 1:4
   for j = 1:7
      c{i,j} = rand(3,2);
   end
end

size(c)
fs = [size(c{1}), size(c)]


%%



%%
clc
cc = c;

cdim = ndims1(c)
ddim = ndims1(data)

csize = size(c)


for cd = 1:cdim
   
   cs = csize(cd+1:end);
   %if length(cs) == 1
   %   cs = [cs,1];
   %end
   cs
   nc = prod(cs)
   
   asgn = num2cell(ones(1, cdim-cd+1));
   asgn{1} = ':';
   asgn
   
   if length(cs) == 1
      ccn = cell(cs,1);
   else
      ccn = cell(cs);
   end
   
   
   for i = 1:nc
      sub = imind2sub(cs, i)
      
      asgn(2:end) = num2cell(sub)
      asgn
      
      size(cc)
      size(cc{1})
      cca = cc(asgn{:})

      
      cat(ddim + cd, cca{:});
      
      ccn{i} = cat(ddim + cd, cca{:});
   end

   cc = ccn;
end

cc = cc{1};

%%

size(cc)
fs


%%
clc
a = 1; b = 1;

c{a,b}
cc(:,:,a,b)





