N = 16;    %number of discretisations
L = 1;       %length of cylinder   
a = 0.01;  %radius of cylinder
m = 3;      %m-point quadrature rule used
eps = 8.854e-12; 

W = L/N;  %discretisation length
disc_l = (W/2:W:1 - W/2);    %discretised vector along length

%rho_l = 2*pi*a*rho_s
% eqn used: 1 = integral( (rho_l * dl) /(4*pi*eps * r), 0, L )

% expressing eqn in Ax = B form
A = ones(N);
B = 4*pi*eps*ones(N,1);
%rho_l: unknown matrix of (N,1)

for i=1:N
    for j=1:N
        A(i,j) = gauss_int( @(yp) inv_dist( 0, W*(2*i - 1)/2, 0,  0, yp, a ), (j-1)*W, j*W, m );
    end
end

rho_l = A \ B;
rho_s = rho_l ./ (2*pi*a);

figure(1);
plot(disc_l, rho_s);
title('Surface Charge distrubution');

figure(2);
imagesc(rho_s'); 
colormap(jet(N));
colorbar;
title('Surface Charge distrubution');

function output = inv_dist(x1,y1,z1, x2,y2,z2)
    dist = sqrt( (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2 );
    output = 1 / dist;
end
