function output = gauss_int(func, lim1, lim2, n)
%Calculates integral value using gaussian quadrature
    %   func: function whose integral is evaluated
    %   lim1,lim2: limits of integration
    %   n: no. points of quadrature evaluation

    mean = (lim1 + lim2) ./ 2;    %scale [lim1,lim2] interval to [-1,1] interval 
    scale = (lim2 - lim1) ./2;
    
    syms x; 
    if n==3
        pts_p = [0, -0.7746, 0.7746];
    else
        pts_p = double(solve(legendreP(n,x) == 0));   %roots of the legendre polynomial
    end
    
    pts = (pts_p .* scale) + mean;   %shift and scale roots
    
    wts = zeros(n);   %generate weights
    for i=1:n
        lagrange_sym = 1;
        for j=1:n
            if j ~= i
                lagrange_sym = conv( lagrange_sym, [1, -pts(j)] ) ./ (pts(i) - pts(j));   %generate lagrange polynomial
            end
        end
        %lagrangeP = matlabFunction(lagrange_sym);
        %wts(i) = integral(lagrangeP, lim1, lim2);
        int_lagrangeP = polyint(lagrange_sym);
        wts(i) = diff( polyval(int_lagrangeP,[lim1, lim2]) );
    end
      
    output = 0;
    for i=1:n
        output = output + func(pts(i)) * wts(i);
    end
    
end
