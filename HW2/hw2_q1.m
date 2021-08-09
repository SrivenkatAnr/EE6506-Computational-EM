fun = @(x) exp(x) ./ sqrt(x.^2);
Nmax = 20;

figure;
a = [0.1:0.01:4.9];
plot(a, fun(a));
title("Plotting the given function")

crct_val = integral(fun,1,4); %The exact value
% Initializing
estimate_val = zeros(4,1);

for i=2:Nmax
    % computing weights and points for evaluation for gauss quadrature
    [pts, wts] = get_wts_pts_for_gauss(1, 4, i);
    % calculating the estimated integral
    estimate_val(i-1) = double( gauss_int(fun,pts, wts) );
end

figure;
plot([2:Nmax], crct_val*ones(Nmax-1,1));
hold on;
plot([2:Nmax], estimate_val);
legend("Exact value (using integral())","Estimate value (using gauss quadrature)",'Location','SouthEast');
xlabel('No. of points of evaluation')
ylabel('Value')
title("Compare Gaussian Quadrature results with Integral results");
