tic
% 1. Defining the parameters
lambda = 0.1;
k0 = 2*pi/lambda;
a = 3*lambda;                   %edge length of pentagon
da = lambda/15;                 %discretization length, please ensure that rem(a,da) = 0
n_sides = 5;                    %number of sides in the object
theta_i_deg = 45;
theta_i = theta_i_deg*(pi/180); %angle of incidence
theta_obj_deg = theta_i_deg+180;%The direction the object is facing
a_ff = 4*lambda;                %far field radius
n_ff = 180;                     %discretizations for far field
tolabs = 1e-9;                  %absolute tolerance in integral
tolrel = 1e-6;                  %relative tolerance in integral
a_sc = 1e-5;                    %radius of deformation (for a segment length of '1')

% 2. Defining the variables 
% 2.1. Related to the object
[test_pt, strt_pt] = get_shape_coords(a, da, n_sides, theta_obj_deg);  

figure;
scatter(test_pt(1,:), test_pt(2,:),'.')
hold on
scatter(strt_pt(1,:), strt_pt(2,:),'.')
axis equal;
legend('delta testing points','starting points for pulse bases')
title('Object boundary')

params = [k0, da, theta_i, a_ff, n_ff, tolabs, tolrel];

% 2.2. Far field points
th_ff = 0:360/n_ff:(360-360/n_ff);
ff_pt = a_ff*[cosd(th_ff); sind(th_ff)];

% 3. Calculating RCS
% For PEC, please set eps_r as Inf
% 3.1. For PEC
eps_r_PEC = Inf;    %The functions compute for PEC when they encounter eps_r = Inf  
fields_bndry = solve_on_boundary(eps_r_PEC, test_pt, strt_pt, params, a_sc); %Calculating the fields on the boundary of the object
RCS_PEC = get_RCS(eps_r_PEC, fields_bndry, ff_pt, strt_pt, params);          %Calculating the RCS

%3.2. For carbon fibre
eps_r_CF = 12-5.5j;
fields_bndry = solve_on_boundary(eps_r_CF, test_pt, strt_pt, params, a_sc);
RCS_CF = get_RCS(eps_r_CF, fields_bndry, ff_pt, strt_pt, params);

toc

%Plotting RCS
figure;
polarplot(th_ff([1:end 1])*pi/180, RCS_PEC([1:end 1]));
title(['RCS for PEC, \theta_i = ',num2str(theta_i_deg),'^\circ'])
figure;
polarplot(th_ff([1:end 1])*pi/180, RCS_CF([1:end 1]));
title(['RCS for carbon fibre, \epsilon_r = ',num2str(eps_r_CF),', \theta_i = ',num2str(theta_i_deg),'^\circ'])


