function output = gauss_int(func, pts, wts)
%Calculates integral value using gaussian quadrature
    %   func: function whose integral is evaluated
    %   pts: points for evaluation
    
    n = length(pts);
    output = 0;
    for i=1:n
        output = output + func(pts(i)) * wts(i);
    end
    
end
