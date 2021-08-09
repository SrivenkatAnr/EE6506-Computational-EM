N = 16;    %number of discretisations
L = 1;       %length of cylinder   
a = 0.01;  %radius of cylinder
m = 3;      %m-point quadrature rule used
R = 10;    %radius of sphere
eps = 8.854e-12; 

W = L/N;  %discretisation length
disc_l = (W/2:W:1 - W/2);    %discretised vector along length

rho_s = 1e-9 .* [0.2698,0.2061,0.1948,0.1882,0.1842,0.1818,0.1803,0.1796,0.1796,0.1803,0.1818,0.1842,0.1882,0.1948,0.2061,0.2698];
rho_l = rho_s .* (2*pi*a);
figure;
plot(disc_l,rho_s);

figure;
n = 6;
theta = pi*(-1:(1/n):1);
phi = (pi/2)*(-1:(1/n):1)';
X = R .* cos(phi)*cos(theta);
Y = R .* cos(phi)*sin(theta);
Z = R .* sin(phi)*ones(size(theta));

V = zeros(2*n + 1);
for i=1:(2*n + 1)
    for j=1:(2*n + 1)
        for k=1:N
            V(i,j) = V(i,j) + rho_s(k) * gauss_int( @(yp) inv_dist( X(i,j), Y(i,j), Z(i,j),  0, yp, 0 ), (k-1)*W, k*W, m );
        end
        
    end
end

surf(X,Y,Z,V);
colormap(jet);
colorbar;
function output = inv_dist(x1,y1,z1,x2,y2,z2)
    dist = sqrt( (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2 );
    output = 1 / dist;
end



