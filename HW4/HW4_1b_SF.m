tic

% Defining the parameters
n_layers = 7;
seq =2;
eps_r = 3.5^2; %relative permitivitty
n = sqrt(eps_r); %refractive index
air_thickness = 1;
ratio = ((sqrt(5) + 1)/2);
k_max = 2*pi;
k_min = 0;
num_pts_per_lyr = 201;

% Discretization lengths for air and medium
DL1 = air_thickness/(num_pts_per_lyr-1);
DL2 = ratio*air_thickness/(num_pts_per_lyr-1);

% The vector with k values
k_vec = k_min:(k_max-k_min)/500:k_max;
%k_vec = 2*pi;
len_vec = length(k_vec);

% Initializing arrays to store ref. and trans. coefficients
tau_arr = zeros(1,len_vec);
ref_arr = zeros(1,len_vec);

% Defining the ref. inde array and the corresponding widths
n_obj_arr = [1 1 1 1 (get_multilayer_eps(seq, n_layers, eps_r)).^0.5 1 1 1 1];
wid_arr = get_width(n_obj_arr, air_thickness, ratio);
wid_arr([1 2 3 4 (end-3) (end-2) (end-1) end]) = air_thickness;

% Combining consecutive objects with the same ref. index (i.e. treating as one, and updating width)
num_objs = length(n_obj_arr);
n_obj_eff_arr = [n_obj_arr(1)];
wid_eff_arr = [wid_arr(1)];
for i = 2:num_objs
    if n_obj_arr(i) == n_obj_arr(i-1)
        wid_eff_arr(end) = wid_eff_arr(end) + wid_arr(i);
    else
        n_obj_eff_arr(end+1) = n_obj_arr(i);
        wid_eff_arr(end+1) = wid_arr(i);
    end
end
num_objs_eff = length(n_obj_eff_arr);

% Effective number of points in each layer, array with position values,
% array with corresponding n values
num_pts_eff_each_lyr = zeros(1,num_objs_eff);
z_arr = 0:DL1:wid_eff_arr(1);
n_arr = ones(1,length(z_arr));
for i = 1:num_objs_eff
    if n_obj_eff_arr(i) == 1
        num_pts_eff_each_lyr(i) = fix(wid_eff_arr(i)/DL1) + 1;
        if i>1
            z_append = sum(wid_eff_arr(1:i-1)):DL1:sum(wid_eff_arr(1:i));
            z_arr = [z_arr, z_append];
            n_arr = [n_arr, ones(1, length(z_append))];
        end
    else
        num_pts_eff_each_lyr(i) = fix(wid_eff_arr(i)/DL2) + 1;
        z_append = sum(wid_eff_arr(1:i-1)):DL2:sum(wid_eff_arr(1:i));
        z_arr = [z_arr, z_append];
        n_arr = [n_arr, n_obj_eff_arr(i)*ones(1, length(z_append))];
    end
end

% Number of equations
num_eq = sum(num_pts_eff_each_lyr);

% Precomputing a few quantities
intN1N1_1 = DL1/3; %integral(N1*N1) for \Delta = DL1
intN1N2_1 = DL1/6;
intN2N2_1 = DL1/3;

intN1N1_2 = DL2/3;
intN1N2_2 = DL2/6;
intN2N2_2 = DL2/3;

% Running the loop over k
for k_id = 1:len_vec
 
        k = k_vec(k_id);
        
        % Initialising
        A = zeros(num_eq);
        b = zeros(num_eq, 1);
        
        % The incident field as a function of position
        Uin = @(z) exp(-1j*k.*z);
        
        % alpha values
        alpha_in = -1j*k;
        alpha_s_l  = -1j*k;
        alpha_s_r = 1j*k;
        
        off_diag_1 = -1/DL1 - k^2*intN1N2_1;
        diag_1 = 2/DL1 - k^2*(intN1N1_1 + intN2N2_1);

        off_diag_2 = -1/DL2 - (k*n)^2*intN1N2_2;
        diag_2 = 2/DL2 - (k*n)^2*(intN1N1_2 + intN2N2_2);
        
        % Defining equations
        id_eq =1;
        A(1,[1 2]) = [(diag_1/2+alpha_s_l) off_diag_1];
        for i = 2:(num_pts_eff_each_lyr(1)-1)
            id_eq = id_eq + 1;
            A(i,[i-1 i i+1]) = [off_diag_1 diag_1 off_diag_1];
        end
        id_eq = id_eq+1;
        A(id_eq, [(id_eq-1) (id_eq) (id_eq+1) (id_eq+2)]) = [off_diag_1 diag_1/2 diag_2/2 off_diag_2];
        b(id_eq) = -alpha_in*Uin(wid_eff_arr(1));
        
        for id_obj = 2:num_objs_eff
            if n_obj_eff_arr(id_obj) ~= 1 % medium
                id_eq = id_eq+1;
                A(id_eq, [id_eq-1 id_eq]) = [-1 1];
                b(id_eq) = Uin(sum(wid_eff_arr(1:id_obj-1)));
                for i = 2:num_pts_eff_each_lyr(id_obj)-1
                    id_eq = id_eq + 1;
                    A(id_eq,[id_eq-1 id_eq id_eq+1]) = [off_diag_2 diag_2 off_diag_2];
                end
                id_eq = id_eq+1;
                if id_obj ~= num_objs_eff
                    A(id_eq, [(id_eq-1) (id_eq) (id_eq+1) (id_eq+2)]) = [off_diag_2 diag_2/2 diag_1/2 off_diag_1];
                    b(id_eq) = alpha_in*Uin(sum(wid_eff_arr(1:id_obj)));
                else % last layer
                    A(id_eq, [(id_eq-1) (id_eq)]) = [off_diag_2 (diag_2/2-alpha_s)];
                end

            else % air
                id_eq = id_eq+1;
                A(id_eq, [id_eq-1 id_eq]) = [1 -1];
                b(id_eq) = Uin(sum(wid_eff_arr(1:id_obj-1)));
                for i = 2:num_pts_eff_each_lyr(id_obj)-1
                    id_eq = id_eq + 1;
                    A(id_eq,[id_eq-1 id_eq id_eq+1]) = [off_diag_1 diag_1 off_diag_1];
                end
                id_eq = id_eq+1;
                if id_obj ~= num_objs_eff
                    A(id_eq, [(id_eq-1) (id_eq) (id_eq+1) (id_eq+2)]) = [off_diag_1 diag_1/2 diag_2/2 off_diag_2];
                    b(id_eq) = -alpha_in*Uin(sum(wid_eff_arr(1:id_obj)));
                else % last layer
                    A(id_eq, [(id_eq-1) (id_eq)]) = [off_diag_1 (diag_1/2-alpha_s_r)];
                end
            end
        end
        
        % Solving the equations
        U = (A\b)';
        
        % Calculating Uin
        Uin_arr = Uin(z_arr);
        
        % Reflection and transmission coefficients
        ref_arr(k_id) = sum(abs( U(1:10) ))/10;
        tau_arr(k_id) = sum(abs( U(end-9:end)+Uin_arr(end-9:end) ))/10;
        
end

% U total in the first and last layers
U_tot_first_lyr = U(1:num_pts_eff_each_lyr(1)) + Uin_arr(1:num_pts_eff_each_lyr(1));
U_tot_last_lyr  = U(end-num_pts_eff_each_lyr(end)+1:end) + Uin_arr(end-num_pts_eff_each_lyr(end)+1:end);
U_tot_bndry_lyrs = [U_tot_first_lyr U_tot_last_lyr];
z_arr_bndry_lyrs = [z_arr(1:num_pts_eff_each_lyr(1)) z_arr(end-num_pts_eff_each_lyr(end)+1:end)];

% Plotting ref. and trans. coefficients
figure;
hold on;
plot(k_vec, tau_arr);
plot(k_vec, ref_arr);
xticks(k_vec(1:fix(len_vec/10):len_vec))
xticklabels(strcat(string(k_vec(1:fix(len_vec/10):len_vec)./pi), '\pi'))
xlabel('k')
hold off;
legend("transmission","reflection");

% Plotting fields
figure
% Magnitude
subplot(2,1,1)
hold on
plot(z_arr, abs(U),'.')
plot(z_arr, abs(Uin_arr))
plot(z_arr_bndry_lyrs, abs(U_tot_bndry_lyrs),'.');
plot(z_arr, n_arr)
legend('abs(U)','abs(Uin)', 'U tot. bndry', 'refr. index')

% Phase
subplot(2,1,2)
hold on
plot(z_arr, (angle(U)),'.');
plot(z_arr, (angle(Uin_arr)));
plot(z_arr_bndry_lyrs, (angle(U_tot_bndry_lyrs)),'.');
plot(z_arr, n_arr)
legend('phase(U)','phase(Uin)', 'U tot. bndry', 'refr. index')

toc