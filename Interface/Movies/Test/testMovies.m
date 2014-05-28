%% test movies


for i = 1:5
   img{i} = rand(25,25);
end

implottiling(img)


%%
images2movie(img, 'test.avi', 'Quality', 50)