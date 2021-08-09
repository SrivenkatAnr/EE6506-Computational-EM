w = 0.05;   %discretisation parameter
x = 2:w:8;  %range for x and y
y = 0:w:6;

[X,Y] = meshgrid(x,y);

Z = exp( -(X-5).^2 - 4.*(Y-3).^2 );   %construct 2d gaussian hill
z_max = max(max(Z));
z_contour = z_max/2;    %find contour at half the height

figure(1);
surf(X,Y,Z);    %plotting the 2d gaussian hill
title("Hill with 2D Gaussian function"); 
hold on;
[~,H] = contour3(X,Y,Z,[z_contour,z_contour]);    %marking contour at half  the height
H.LineWidth = 8;
axis equal;
hold off;

figure(2);
surf(X,Y,Z);    %plotting the 2d gaussian hill
title("Normals at the contour");    
hold on;
[C,H] = contour3(X,Y,Z,[z_contour,z_contour]);    %plot contour
H.LineWidth = 8;

P = C(:,2:101)';
Eq_P = curvspace(P,25)';    %find equispaced points on the contour

[nx,ny,nz] = find_normal(Z,w);    %defined function to calculate normal
filter = filterxy(X,Y,Eq_P);
quiver3(X, Y, Z, nx .* filter, ny .* filter, nz .* filter);    %plotting normals at equispaced points
axis equal;
hold off;

%{ 
% Validation for finding the normal done using the inbuilt surfnorm function
figure(3);
surf(X,Y,Z);
title("Validation: Normals found from surfnom"); 
hold on;
[~,C] = contour3(X,Y,Z,[z_contour,z_contour]);
C.LineWidth = 8;
[U,V,W] = surfnorm(X,Y,Z);
quiver3(X,Y,Z,U,V,W);
axis equal;
hold off;
%}

function filter = filterxy(X,Y,q)    %to overlay equispaced points on the predefined mesh
    tol = 0.026;
    filter = zeros(size(X));
    %calculate nearest points on the mesh for each one of the equispaced
    %points. then do elementwise OR for individual filters to get net filter
    for i=1:size(q')
        temp_filter = X < (q(1,i) + tol) & X > (q(1,i) - tol) & Y < (q(2,i) + tol) & Y > (q(2,i) - tol);
        filter = filter | temp_filter;
    end    
end

function [nx,ny,nz] =  find_normal(g,w)    %to find normal for a surface g
    %g = z - f(x,y),  normal(g) = 1 - grad(f)
    [nx,ny] =  gradient(g .* (-1/w));    
    nz = ones([1 + (6/w), 1 + (6/w)]);
end
