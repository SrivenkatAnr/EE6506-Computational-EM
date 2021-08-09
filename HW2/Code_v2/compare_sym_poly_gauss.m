fun = @(x) exp(x) ./ sqrt(x.^2);
Nmax = 25;

crct_val = integral(fun,1,4); %The exact value
% Initializing
val1 = zeros(Nmax-1,1);
val2 = zeros(Nmax-1,1);
val3 = zeros(Nmax-1,1);

syms x;
for n=2:Nmax
    % computing weights and points for evaluation for gauss quadrature
    [pts, wts_sym] = get_wts_pts_for_gauss(1, 4, n);    
    wts_pol = zeros(n,1);   
    wts_sympol = zeros(n,1);
    
    for i = 1:n
        lagrange_sym = 1;
        lagrange_pol = 1;
        for j = 1:n
            if i ~= j
                lagrange_sym = lagrange_sym*(x - pts(j)) / (pts(i) - pts(j));
                lagrange_pol = conv( lagrange_pol, [1, -pts(j)] ) ./ (pts(i) - pts(j)); 
            end
        end
            lagrange_sympol = sym2poly(lagrange_sym);
            
            int_lagrange_pol = polyint(lagrange_pol);
            wts_pol(i) = diff( polyval(int_lagrange_pol,[1, 4]) );
            
            int_lagrange_sympol = polyint(lagrange_sympol);
            wts_sympol(i) = diff( polyval(int_lagrange_sympol,[1, 4]) );            
    end    
        
    % calculating the estimated integral
    val1(n-1) = double( gauss_int(fun,pts, wts_sym) );
    val2(n-1) = double( gauss_int(fun,pts, wts_pol) );
    val3(n-1) = double( gauss_int(fun,pts, wts_sympol) );
end

figure;
plot([2:Nmax], crct_val*ones(Nmax-1,1));
hold on;
plot([2:Nmax], val1,'o','markersize',10);
plot([2:Nmax], val2,'x','markersize',12);
plot([2:Nmax], val3,'.','markersize',12);
legend("Exact value (using integral())","Full symbolic","Full Polynomial","Symbolic+Polynomial",'Location','NorthWest');
xlabel('No. of points of evaluation')
ylabel('Value')
title("Compare Different Gaussian Quadrature results with Integral results");
