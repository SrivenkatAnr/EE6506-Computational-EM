fun = @(x) exp(x) ./ sqrt(x.^2);
Nmax = 20;

figure(1);
a = [0.1:0.01:4.9];
plot(a, fun(a));
title("Plotting the given function")
% A polynomial approx of order 3 or more should decently fit by intuition. 

crct_val = integral(fun,1,4)
estimate_val = zeros(4,1);
for i=2:Nmax
    estimate_val(i-1) = double( gauss_int(fun,1,4,i) );
end

figure(2);
plot([2:Nmax], crct_val*ones(Nmax-1,1));
hold on;
plot([2:Nmax], estimate_val);
legend("Estimated value","Accurate value",'Location','SouthEast');
title("Compare Gaussian Quadrature results with Integral results");