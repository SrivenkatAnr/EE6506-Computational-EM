function [pts, wts] = get_wts_pts_for_gauss(lim1, lim2, n)
    %   lim1,lim2: limits of integration
    %   n: no. points of quadrature evaluation
    mean = (lim1 + lim2) ./ 2;    %scale [lim1,lim2] interval to [-1,1] interval 
    scale = (lim2 - lim1) ./2;
    
    syms x; 
    pts_p = double(solve(legendreP(n,x) == 0));   %roots of the legendre polynomial
    
    pts = (pts_p .* scale) + mean;   %shift and scale roots
    
    wts = zeros(n,1);   %generate weights
    for i = 1:n
        lagrange_sym = 1;
        for j = 1:n
            if j ~= i
                lagrange_sym = lagrange_sym*(x - pts(j))/(pts(i) - pts(j));    %Construct lagrange polynomial
            end
        end        
        lagrange_pol = coeffs(lagrange_sym,'all');    %get the coefficients of lagrangeP with high precision
        int_lagrange_pol = poly2sym(polyint(lagrange_pol));    %integrate the lagrangeP and convert into symbolic form
        wts(i) = subs(int_lagrange_pol,lim2) - subs(int_lagrange_pol,lim1);     %calc weights
    end

end
