function L_coeff = scalar_basis(pt_arr)
        %L_{i}(x) is linear fn of (xi,yi), i.e.,  (ax + by + c)/k 
        %this returns [a,b,c,k]
        
        a = pt_arr(2,2) - pt_arr(3,2);
        b = pt_arr(3,1) - pt_arr(2,1);
        c = ( pt_arr(2,1) * pt_arr(3,2) ) - ( pt_arr(3,1) * pt_arr(2,2) );
        k = a*pt_arr(1,1) + b*pt_arr(1,2) + c;
        
        L_coeff = [a, b, c, k];     %express basis fn as array of coeff