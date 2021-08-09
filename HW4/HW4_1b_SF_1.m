tic

n_layers = 1;
seq =1;
eps_r = 3.5^2; %relative permitivitty
n = sqrt(eps_r); %refractive index
air_thickness = 1;
ratio = ((sqrt(5) + 1)/2);
k_max = 2*pi;
k_min = 0;
num_pts_per_lyr = 201;
DL1 = air_thickness/(num_pts_per_lyr-1);
DL2 = ratio*air_thickness/(num_pts_per_lyr-1);

k_vec = k_min:(k_max-k_min)/500:k_max;
%k_vec = 2*pi;
len_vec = length(k_vec);
tau_arr = zeros(1,len_vec);
ref_arr = zeros(1,len_vec);

n_obj_arr = [1 1 (get_multilayer_eps(seq, n_layers, eps_r)).^0.5 1 1];
wid_arr = get_width(n_obj_arr, air_thickness, ratio);
%wid_arr([1 2 (end-1) end]) = air_thickness;

% n_obj_arr = [1 1.5];
% wid_arr = [5 15];

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

num_eq = sum(num_pts_eff_each_lyr);

for k_id = 1:len_vec
 
        k = k_vec(k_id);
        
        A = zeros(num_eq);
        b = zeros(num_eq, 1);

        Uin = @(z) exp(-1j*k.*z);

        alpha_in = -1j*k;
        alpha_s  = 1j*k;

        intN1N1_1 = DL1/3;
        intN1N2_1 = DL1/6;
        intN2N2_1 = DL1/3;

        intN1N1_2 = DL2/3;
        intN1N2_2 = DL2/6;
        intN2N2_2 = DL2/3;

        off_diag_1 = -1/DL1 - k^2*intN1N2_1;
        diag_1 = 2/DL1 - k^2*(intN1N1_1 + intN2N2_1);

        off_diag_2 = -1/DL2 - (k*n)^2*intN1N2_2;
        diag_2 = 2/DL2 - (k*n)^2*(intN1N1_2 + intN2N2_2);
        
        %left most interface
        id_eq =1;
        A(1,[1 2]) = [(diag_1/2+alpha_s) off_diag_1];
        for i = 2:(num_pts_eff_each_lyr(1))
            id_eq = id_eq + 1;
            A(i,[i-1 i i+1]) = [off_diag_1 diag_1 off_diag_1];
        end

        for id_obj = 2:(num_objs_eff-1)
            if n_obj_eff_arr(id_obj) ~= 1    %object made of material N
                for i = 1:num_pts_eff_each_lyr(id_obj)
                    id_eq = id_eq + 1;
                    A(id_eq,[id_eq-1 id_eq id_eq+1]) = [off_diag_2 diag_2 off_diag_2];
                end

            else    %object made of air
                for i = 1:num_pts_eff_each_lyr(id_obj)
                    id_eq = id_eq + 1;
                    A(id_eq,[id_eq-1 id_eq id_eq+1]) = [off_diag_1 diag_1 off_diag_1];
                end
            end
        end
         
        %right most interface
        for i = 1:(num_pts_eff_each_lyr(end)-1)
            id_eq = id_eq + 1;
            A(id_eq,[id_eq-1 id_eq id_eq+1]) = [off_diag_1 diag_1 off_diag_1];
        end
        id_eq = id_eq + 1;
        A(id_eq, [(id_eq-1) (id_eq)]) = [off_diag_1 (diag_1/2-alpha_s)];

        U = (A\(b))';

        Uin_arr = Uin(z_arr);
        
        ref_arr(k_id) = sum(abs( U(1:10) ))/10;
        tau_arr(k_id) = sum(abs( U(end-9:end)+Uin_arr(end-9:end) ))/10;
        
end

figure;
hold on;
plot(k_vec, tau_arr);
plot(k_vec, ref_arr);
xticks(k_vec(1:fix(len_vec/10):len_vec))
xticklabels(strcat(string(k_vec(1:fix(len_vec/10):len_vec)./pi), '\pi'))
xlabel('k')
hold off;
legend("transmission","reflection");

figure
subplot(2,1,1)
hold on
plot(z_arr, abs(U),'.')
plot(z_arr, abs(Uin_arr))
plot(z_arr, n_arr)
legend('abs(U)','abs(Uin)', 'refr. index')

subplot(2,1,2)
hold on
plot(z_arr, unwrap(angle(U)),'.');
plot(z_arr, unwrap(angle(Uin_arr)));
plot(z_arr, n_arr)
legend('phase(U)','phase(Uin)', 'refr. index')

toc
