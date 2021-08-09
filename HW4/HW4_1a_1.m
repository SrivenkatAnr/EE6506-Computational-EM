tic

n = 7;    %number of layers
seq =1;  %choose structure 1 or 2
eps_n = 3.5^2; %relative permitivitty
air_thickness = 1;
ratio = ((sqrt(5) + 1)/2);    %ratio of thicknesses of material layer and air gap
k_vec = 0:2*pi/1000:2*pi;  %sweep k instead of lambda to check for periodicity

len_vec = length(k_vec);
tou_arr = zeros(1,len_vec);
ref_arr = zeros(1,len_vec);

eps_arr = get_multilayer_eps(seq, n, eps_n);   
wid_arr = get_width(eps_arr, air_thickness, ratio);
len = size(eps_arr');

tou = @(eps1, eps2) 2*sqrt(eps1) ./ (sqrt(eps1) + sqrt(eps2));     %transmission coeff. from fresnel eqn
ref = @(eps1, eps2)  (sqrt(eps1) - sqrt(eps2)) ./ (sqrt(eps1) + sqrt(eps2));      %reflection coeff. from fresnel eqn
 
I_mat  =@(t,r) (1/t) .* [1 r; r 1;];    %imposes constraints at every interface 
delta = @(eps,wid, k) k*sqrt(eps)*wid;
P_mat = @(delta) [exp(1i*delta) 0; 0 exp(-1i*delta)];    %includes phase difference across a layer

eps_arr(end+1) = 1;

for k_id = 1:len_vec
        %get net transfer matrix
        T_mat = I_mat( tou(1, eps_arr(1)), ref(1, eps_arr(1)) );    
        for i = 1:len
              T_mat = T_mat * P_mat( delta(eps_arr(i), wid_arr(i), k_vec(k_id)) ) * I_mat( tou(eps_arr(i), eps_arr(i+1)), ref(eps_arr(i), eps_arr(i+1)) );     
        end
        
        %get tou, ref from T matrix 
        net_tou = 1 / abs(T_mat(1,1));
        net_ref = abs( T_mat(2,1) / T_mat(1,1) );

        tou_arr(k_id) = net_tou;
        ref_arr(k_id) = net_ref;

end

toc

figure;
hold on;
plot(k_vec, tou_arr);
plot(k_vec, ref_arr);
xticks(k_vec(1:fix(len_vec/10):len_vec))
xticklabels(strcat(string(k_vec(1:fix(len_vec/10):len_vec)./pi), '\pi'))
xlabel('k')
hold off;
legend("transmission","reflection");
