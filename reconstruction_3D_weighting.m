function [vol_r] = reconstruction_3D_weighting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3D reconstraction of the original volume using its 2d projections

% written by Hassan keshvari Khojasteh...

%%%%%%%%%%%%%%%%%%%% load the projections and prepare them %%%%%%%%%%%%%%%%

names = dir(fullfile('projecs/*.mat'));
vol_rec = [];

cnt = 1;
for ii = 1:length(names)

    name_ = names(ii).name;
    name_ = append('projecs/', name_);
    
    img_proj = load(name_).m_ref;
    [d_1, d_2] = size(img_proj);

    ss_sum = sum(img_proj);
    [~, l] = find(ss_sum == 0);
    
    if isempty(l)
        shiftx = 0;
    else
        shiftx = l(length(l));
    end    
    
    img = zeros(size(img_proj));
    img(:, 1:d_2-shiftx) = img_proj(:, shiftx+1:d_2);
    theta = str2double(name_(16:length(name_)-4));
   
    %%%%%%%%%%%%%%%%%%%% reconstruct the volume %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isempty(vol_rec)
        nz = 32; 
        nx = d_1;
        ny = d_2;
        
        vol_rec = zeros(d_1, d_2, nz);        
        midd = (length(vol_rec)+1)/2;

        xx = linspace(1, d_1, d_1);
        yy = linspace(1, d_2, d_2);
        zz = linspace(1, nz, nz);
        xx = xx - midd;
        yy = yy - midd;
        zz = zz - midd;
    end   
    
    M = rot_mat(theta, 0 , 0);
      
    for iz = 1:nz
        for iy= 1:ny
            for ix= 1:nx
                rr = M*[xx(ix); yy(iy); zz(iz)] + midd;
                val = interp2(img, rr(1), rr(2), 'linear');
                if isnan(val)
                   val = 0;
                end   
                vol_rec(xx(ix) +midd , yy(iy)+ midd, zz(iz)+ midd) = vol_rec(xx(ix) +midd , yy(iy)+ midd, zz(iz)+ midd) + val;
            end
        end    
    end
    
    sprintf('data reconstructed using %d projections', cnt)
    cnt = cnt + 1;
end


figure(1)
isosurface(vol_rec);
title('reconstructed data before weighting')
xlabel('x')
ylabel('y')
zlabel('z')  

%%%%%%%%%%%%%%% weighting the reconstructed volume %%%%%%%%%%%%%%%%%%%%%%%%

fvol = fftshift(fft(fft(fft(vol_rec,[],1),[],2),[],3));

[x,y,z]=ndgrid(1:d_1, 1:d_2, 1:nz);
midd = (length(vol_rec)+1)/2;
x = reshape(x, d_1*d_2*nz, 1) - midd;
y = reshape(y, d_1*d_2*nz, 1) - midd;
z = reshape(z, d_1*d_2*nz, 1) - midd;
t = [x, y, z];

angls = round(linspace(-90, 90, length(names)), 4);
two_a = sqrt(x(d_1).^2*3);
N = length(angls);

w = (zeros(size(vol_rec)));

for j = 1:N
    rr = rot_mat(angls(j), 0, 0)*t';
    z_j = rr(3, :) + midd;
    w = w + reshape(two_a*sinc(two_a*pi*z_j), [d_1, d_2, nz]);
end

w=1./(1+w);
min_w=min(min(min(w)));
max_w=max(max(max(w)));
w=(w-min_w)./(max_w-min_w);

threshold = mean(mean(mean(w)));
ind = find(abs(w)<threshold); 
w(ind) = threshold;   

radius = 1.5;
filter = sqrt(reshape(x, d_1, d_2, nz).^2 + reshape(y, d_1, d_2, nz).^2 + reshape(z, d_1, d_2, nz).^2) <= radius;
w=w.*filter;
    
out = fvol.*w;

%%%%%%%%%%%%%%%%%%%% back to the real space %%%%%%%%%%%%%%%%%%%%%%%%%%%

vol_r = real(ifft(ifft(ifft(ifftshift(out),[],1),[],2),[],3)); 

figure(2)
isosurface(vol_r);
title('reconstructed data (weighted)')
xlabel('x')
ylabel('y')
zlabel('z')  

r = load ('volume');
r = r.r;

figure(3)
subplot(1,2, 1)
nz_s = randi(length(r), 1);
slice = r(:, :, nz_s);
slice_m = mean(mean(slice));
slice_std = std(std(slice));
slice = (slice - slice_m)./slice_std; 
imshow(slice)
title('sample slice of orig vol')

subplot(1,2, 2)
slice = vol_r(:, :, nz_s);
slice_m = mean(mean(slice));
slice_std = std(std(slice));
slice = (slice - slice_m)./slice_std; 
imshow(slice)
title('sample slice of recons vol')
end


