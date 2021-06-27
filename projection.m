function [] = projection 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generating the dataset of 2d projections

% written by Hassan keshvari Khojasteh...

%%%%%%%%%%%%%%%% Let's plot the 3d structure/2d slice of molecule %%%%%%%%%

r = load ('volume');
r = r.r;
[d_1, d_2, d_3] =size(r);

figure(1)
isosurface(r);
title('original data')
xlabel('x')
ylabel('y')
zlabel('z')       
hold on

figure(2)
nz = randi(length(r), 1);
slice = r(:, :, nz);
slice_m = mean(mean(slice));
slice_std = std(std(slice));
slice = (slice - slice_m)./slice_std; 
imshow(slice)
title('sample slice')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% padding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_pad = zeros(size(r)*3);
r_pad(d_1+1:2*d_1, d_2+1:2*d_2, d_3+1:2*d_3) = r;

%%%%%%%%%%%%%%%%%%%%% rotating the volume. e.g. Cryo-em %%%%%%%%%%%%%%%%%%%

[dr_1, dr_2, dr_3] =size(r_pad);
[x, y, z] = meshgrid(1:dr_1, 1:dr_2, 1:dr_3);

x = reshape(x, dr_1*dr_2*dr_3, 1);
y = reshape(y, dr_1*dr_2*dr_3, 1);
z = reshape(z, dr_1*dr_2*dr_3, 1);
midd = (length(r_pad)+1)/2;

t = [x, y, z];
t = t - midd;
cnt = 1;

for th = round(linspace(-90, 90, 100), 4)
    
    % Starting the roation using the rotation matrix. The angles should be in
    % degree
    theta = th;
    phi = 0;
    psi = 0; 

    M = rot_mat(theta, phi, psi);
    rot = M\t' + midd;

    x_ref = rot(1,:);
    y_ref = rot(2,:);
    z_ref = rot(3,:);

    r_rotated = interp3(r_pad, x_ref, y_ref, z_ref, 'linear');
    r_rotated = reshape(r_rotated, size(r_pad));
    
    idx=find(abs(r_rotated)>0);
    [mx, my, mz] = ind2sub(size(r_rotated), idx);
    r_rotated = r_rotated(min(mx):max(mx), min(my):max(my), min(mz):max(mz));
    
    %%%%%%%%%%%%%%%%%%% 2D projections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % xy plane projection
    d_1w = max(mx) - min(mx) + 1;
    d_2w = max(my) - min(my) + 1;
    d_3w = max(mz) - min(mz) + 1;

    m = zeros(d_1w, d_2w);
    for i = 1:d_3w
         m = m + r_rotated(:, :, i);
    end
    
    m = imresize(m, [d_1, d_2]);
    m_noise = awgn(m, -8,'measured');% adding the noise
    m_scale = m_noise * 1e-3;

    shift = randi(d_2/2, 1);
    m_scale_sh = circshift(m_scale, shift, 2);
    m_ref = zeros(size(m_scale_sh));
    m_ref(:, shift+1:d_2) = m_scale_sh(:, shift+1:d_2);
    
    figure(3)
    imshow(m_ref);
    title('projected image')
    hold on
    
    if exist('projecs', 'file') == 0
        mkdir 'projecs'
    end   
    
    sprintf('%d images generated!', cnt)
    
%     save(sprintf('projecs/theta_m = %+2.4f.mat', th), 'm') 
%     save(sprintf('projecs/theta_scale = %+2.4f.mat', th), 'm_scale') 
    save(sprintf('projecs/theta = %+2.4f.mat', th), 'm_ref')
    
    cnt = cnt + 1;
end

end




