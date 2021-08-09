Nc = 100;   %number of discretisations of cylinder
Ns = 100;   %number of discretisations of sphere
L = 1;       %length of cylinder   
a = 0.01;  %radius of cylinder
R = 10;    %radius of sphere
eps = 8.854e-12; 
m = 3;      %m-point quadrature rule used

W = L/Nc;  %discretisation length
disc_l = (W/2:W:1 - W/2);    %discretised vector along length

%Calculating points of evaluation and weights for 'm' point gauss
%quadrature
[pts, wts] = get_wts_pts_for_gauss(0, W, m);

%rho_l = 2*pi*a*rho_s
% eqn used: 1 = integral( (rho_l * dl) /(4*pi*eps * r), 0, L )

% expressing eqn in Ax = B form
A = ones(Nc);
B = 4*pi*eps*ones(Nc,1);
%rho_l: unknown matrix of (Nc,1)

for i=1:Nc
    for j=1:Nc
        A(i,j) = gauss_int( @(yp) inv_dist( 0, W*(2*i - 1)/2, 0,  0, yp, a ), pts + (j-1)*W, wts );
    end
end

rho_l = A \ B; %solving for rho_l
rho_s = rho_l ./ (2*pi*a);

%Plotting
figure;
plot(disc_l, rho_s);
xlabel('y')
ylabel('Surface charge density (C/m^2)')
title(['Surface Charge distrubution, N_c = ',num2str(Nc)]);

figure;
imagesc(rho_s'); 
colormap(jet(Nc));
colorbar;
title(['Surface Charge distrubution, N_c = ', num2str(Nc)]);

% Defining points on sphere
phi = pi * (-1:(1/Ns):1);
theta = (pi/2) * (-1:(1/Ns):1)';
X = R .* cos(theta) * cos(phi);
Y = R .* cos(theta) * sin(phi);
Z = R .* sin(theta) * ones(size(phi));

% Matrix for the potential on sphere
V = zeros(2*Ns + 1);

% Calculating the potential over the sphere
for i=1:(2*Ns + 1)
    for j=1:(2*Ns + 1)
        for k=1:Nc
            V(i,j) = V(i,j) + rho_l(k) * gauss_int( @(yp) inv_dist( X(i,j), Y(i,j), Z(i,j),  0, yp, 0 ), pts + ((k-1)-Nc/2)*W, wts );
        end        
    end
end

% Plotting
figure;
colormap(jet);
surf(X,Y,Z,V, 'EdgeColor','interp');
colorbar;
axis square;
xlabel('x')
ylabel('y')
zlabel('z')
title("Potential distribution on sphere of R=10");

function output = inv_dist(x1,y1,z1, x2,y2,z2)
    dist = sqrt( (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2 );
    output = 1 / dist;
end
