Nc = 20;   %number of discretisations of cylinder
Ns = 18;   %number of discretisations of sphere
L = 1;       %length of cylinder   
a = 0.01;  %radius of cylinder
R = 10;    %radius of sphere
m = 3;      %m-point quadrature rule used
eps = 8.854e-12; 

W = L/Nc;  %discretisation length
disc_l = (W/2:W:1 - W/2);    %discretised vector along length

%rho_l = 2*pi*a*rho_s
% eqn used: 1 = integral( (rho_l * dl) /(4*pi*eps * r), 0, L )

% expressing eqn in Ax = B form
A = ones(Nc);
B = 4*pi*eps*ones(Nc,1);
%rho_l: unknown matrix of (Nc,1)

for i=1:Nc
    for j=1:Nc
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
colormap(jet(Nc));
colorbar;
title('Surface Charge distrubution');

phi = pi * (-1:(1/Ns):1);
theta = (pi/2) * (-1:(1/Ns):1)';
X = R .* cos(theta) * cos(phi);
Y = R .* cos(theta) * sin(phi);
Z = R .* sin(theta) * ones(size(phi));

V = zeros(2*Ns + 1);

for i=1:(2*Ns + 1)
    for j=1:(2*Ns + 1)
        for k=1:Nc
            V(i,j) = V(i,j) + rho_l(k) * gauss_int( @(yp) inv_dist( X(i,j), Y(i,j), Z(i,j),  0, yp, 0 ), (k-1)*W, k*W, m );
        end        
    end
end

figure(3);
colormap(jet);
surf(X,Y,Z,V);
colorbar;
axis square;
title("Potential distribution on sphere of R=10");

function output = inv_dist(x1,y1,z1, x2,y2,z2)
    dist = sqrt( (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2 );
    output = 1 / dist;
end
