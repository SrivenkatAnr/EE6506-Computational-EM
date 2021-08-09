%set chosen points for arbitrary triangle
pt1 = [0,0]; 
pt2 = [1,0]; 
pt3 = [0,1];
pt_arr = [pt1; pt2; pt3];

%define scalar basis fn coefficient arrays
L1 = scalar_basis(pt_arr);
L2 = scalar_basis( circshift(pt_arr,2) );
L3 = scalar_basis( circshift(pt_arr,1) );

L_xy = @(x, y, L) (L(1).*x + L(2).*y + L(3)) ./ L(4); %for testing scalar basis functions
grad_L = @(L) [L(1)./L(4), L(2)./L(4), 0];    %to get gradient of L

grad_L1 = grad_L(L1);
grad_L2 = grad_L(L2);
grad_L3 = grad_L(L3);

L_arr = [L1; L2; L3];
grad_L_arr = [grad_L1; grad_L2; grad_L3];

%convention followed in transformation: pt1 -> (0,0), pt2 -> (1,0), pt3 -> (0,1)
%transformation: x = x1 + (x2-x1)u + (x3-x1)v, y = y1 + (y2-y1)u + (y3-y1)v
% (x2-x1) = L3(2), (x3-x1) = -L2(2), (y2-y1) = -L3(1), (y3-y1) = L2(1)

J = [pt2(1) - pt1(1), pt3(1) - pt1(1); pt2(2) - pt1(2), pt3(2) - pt1(2)];    %Jacobian matrix

%coordinate transformation from (x,y) to (u,v)
x_uv = @(u, v) pt_arr(1,1) + ( pt_arr(2,1) - pt_arr(1,1) ) .* u +  ( pt_arr(3,1) - pt_arr(1,1) ) .* v;    
y_uv = @(u, v) pt_arr(1,2) + ( pt_arr(2,2) - pt_arr(1,2) ) .* u +  ( pt_arr(3,2) - pt_arr(1,2) ) .* v;
L_uv = @(u,v,L) L_xy( x_uv(u,v), y_uv(u,v), L );

A = zeros(3);
B = zeros(3);

y_max = @(x) 1 - x;
for i=1:3
        for j=1:3                
                A(i,j) = integral2( @(u,v) ( L_uv(u,v, L_arr(i,:)) .* L_uv(u,v, L_arr(j,:)) ), 0,1, 0,y_max);
                B(i,j) = dot( grad_L_arr(i,:), grad_L_arr(j,:) ) .* 0.5;
        end
end

A = A .* det(J)  
B = B .* det(J)

%validating the basis functions
x = [0:0.1:6];
y = [0:0.1:6];
[X,Y] = meshgrid(x,y);
Z1 = L_xy(X,Y,L1);    
surf(X,Y,Z1);    %plot L1
hold on;
Z2 = L_xy(X,Y,L2);
surf(X,Y,Z2);    %plot L2
Z3 = L_xy(X,Y,L3);
surf(X,Y,Z3);    %plot L3
xlim([0 max(pt_arr(:,1))]);
ylim([0 max(pt_arr(:,2))]);
zlim([0 1]);

