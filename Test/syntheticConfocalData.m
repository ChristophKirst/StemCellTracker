function [img, ctr] = syntheticConfocalData()
%
% [img, ctr] = syntheticConfocalData()
%
% Create a random array of non overlapping spheres, subject each to an affine
% deformation with random directions but fixed stretching. The ellipses
% resulting from the deformation may overlap. Intensity within each sphere is a
% defined function of radius. Spheres may extend out of the imaged box in x,y
% but not in z. Routine in_sphere() defines intensity in sphere
%
% img = density image in 3D returned. 
% ctr = [ncells, 3] array of centers
%   Note internal array lbl which is 3D label matrix of which cell. Max(lbl)
% == ncells+1 and denotes all overlaps of ellipse. Note graphics at end.
%

rad0 = 2; % radius of sphere, grid units integer+>0
dmin = 20;  % min allowed center to center distance of spheres
% dmin = max(dmin, 2*rad0);
% slab = [256, 256, 48];     % 3D mesh for data

slab = [160, 180, 40];     % 3D mesh for data

slab_r = [slab(1:2), slab(3) - 2*rad0];  % slab dims for centers of cells
volfrac = 0.50;     % approximate vol fraction for spheres, may not be reachable with random addition of pts
ncell0 = 400; % round( volfrac*prod(slab)/(pi*dmin^3/6) );
scl = [1.3,1, 0.8];   % affine scale change applied to r, vol preserving
scl = scl/prod(scl)^(1/3);

% allocate space
img  = single(zeros(slab));     % the synthetic image returned
occ  = false(slab);     % to keep track of excluded volume, could have used false(slab_r)
lbl  = uint16(zeros(slab));  % color labels for cells
lbl_mx = 2^16 - 1;  % special label for overlaps of cells -> white
ctr  = zeros(ncell0, 3);    % centers of spheres.
r_dmin = in_sphere(dmin);    % all pts in sphere radius<dmin centered at 0.
[r_rad0, density] = in_sphere(rad0);

nc = 1;
for i = 1:(10*ncell0) % typically  
    if i < 2*ncell0  % random addition with test at the beginning
        r0 = ceil(rand(1,3) .* slab_r) + [0,0,rad0];
        if occ(r0(1),r0(2),r0(3))
            continue
        end
    else    % when only a few allowed spots find them and sample
        indx = find(~occ);  lindx = length(indx);
        % no more space left when slab_b grid all occupied
        if lindx < 2   
            break;
        end
        
        
        indx1 = randsample(lindx, min(lindx, 20));
        indx1 = indx(indx1);

        for iter = indx1'
            [r0(1), r0(2), r0(3)] = ind2sub(slab, iter);
            % nuclei are confined to the slab in z direction.
            if rad0 <= r0(3) && r0(3) <= slab(3)-rad0+1
                break;
            end
        end
    end
    ctr(nc, :) = r0;
    % update matrix of excluded sites
    all_r = limit_pts(r_dmin + ones(size(r_dmin,1),1)*r0, slab);
    indx = sub2ind(slab, all_r(:,1), all_r(:,2), all_r(:,3) );
    occ(indx) = 1;
    % affine distort the sphere and add density to simulated image.
    all_r = limit_pts( rand_strain(scl, r_rad0) + ones(size(r_rad0,1),1)*r0, slab);
    indx = sub2ind(slab, all_r(:,1), all_r(:,2), all_r(:,3) );
    % can randomize level of signal from cell here, by multiplying density
    facc = 10*rand();
    img(indx) = max(img(indx), facc*density);
    % update label matrix
    overlap = lbl(indx) > 0;
    lbl(indx(overlap)) = lbl_mx;
    lbl(indx(~overlap)) = nc;
    nc = nc + 1;
    if nc > ncell0
        break
    end
end

free_space = sum(sum(sum(~occ(:,:,rad0:(slab(3)-rad0+1)) )));
fprintf('added %d cells, radius= %d to slab= %d %d %d, vol frac= %5.3f, free pts for cells= %d\n',...
    nc-1, rad0, slab, nc*4*pi*rad0^3/3/prod(slab), free_space );

% clean up output
ctr(nc:end,:) = [];
disk = strel('disk', 1);
if rad0 > 5
 %   disk = strel('disk',2); causes more blob mergers
end
for k = 1:slab(3)
    img(:,:,k) = imclose(img(:,:,k), disk);
    img(:,:,k) = imfill(img(:,:,k), 'holes');  % not needed if strel radius=2
    lbl(:,:,k) = imclose(lbl(:,:,k), disk);
end

% color map for display.
% lbl(lbl==lbl_mx) = nc;
% cmap = colorcube(nc-2);
% cmap = [0 0 0; cmap; 1 1 1];    % no cell-> black, overlap -> white
% k = 25;
% rgb = ind2rgb(lbl(:,:,k), cmap);
% figure, imshow(rgb,[]);
% title('color map at arb z-slice for test');
return
%%%%%%% end of main

function rr = rand_strain(scl, rr)
% n? are row vectors and become image of x,y,z components of rr input.
n1 = rand_dir();   
n2 = rand_dir();
n2 = n2 - dot(n1,n2)*n1;
n2 = n2/norm(n2,2);
n3 = cross(n1,n2);
rr = rr(:,1)*scl(1)*n1 + rr(:,2)*scl(2)*n2 + rr(:,3)*scl(3)*n3;
rr = round(rr);

function nn = rand_dir()
% random unit vector
rr = ones(1,3);
while rr * rr' > 1
    rr = 2*rand(1,3) - 1;
end
nn = rr/norm(rr, 2);

function [rr, density] = in_sphere(rad0)
% return all lattice coordinates of points dst < rad0 of origin, along with
% density = what gets added to img to simulate DAPI
npt = 2*rad0;
rr = zeros(npt^3,3);
density = zeros(npt^3,1);
ptr = 1;
for i = 0:npt
    for j = 0:npt
        for k = 0:npt
            rad2 = (i-rad0)^2 + (j-rad0)^2 + (k-rad0)^2;
            if rad2 >= rad0^2
                continue
            else
                rr(ptr, :) = [i,j,k] - rad0;
                density(ptr) = sqrt(rad0^2 - rad2);
                ptr = ptr + 1;
            end
        end
    end
end
rr(ptr:end,:) = [];
density(ptr:end) = [];

function rr = limit_pts(rr, slab)
npts = size(rr,1);
rr = max(rr, 1);
rr = min(rr, ones(npts,1)*slab );

