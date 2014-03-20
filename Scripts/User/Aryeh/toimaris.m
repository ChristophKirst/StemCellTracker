%% load surfaces and transfer to imaris


%%
load('Z:\140305_RUES2_36hBMP4_Bra_Snail_Sox2_surfaces_imaris.mat')

%%
surf = surfaces{1};
fac = surfaces{2};
norm = surfaces{3};


%%
   
nset = 10;
%nset = size(surf)
sfset = surf(1:nset);
fcset = fac(1:nset);
nmset = norm(1:nset);

   
%%  
imarissetsurface('Aryeh Segments', sfset, fcset, nmset);
