function [M] = rot_mat(theta, phi, psi)
% Rotation around the x axis
Rx = [1      0            0     ;
      0 cosd(theta) -sind(theta);
      0 sind(theta)  cosd(theta)];

% Rotation around the y axis
Ry = [ cosd(phi)  0 sind(phi); 
          0       1     0    ;
      -sind(phi)  0 cosd(phi)];


% Rotation around the z axis
Rz = [cosd(psi)   -sind(psi)   0;
      sind(psi)    cosd(psi)   0; 
           0              0    1];


M = Rx*Ry*Rz;
end
