function RCS = get_RCS(eps_r, fields_bndry, ff_pt, strt_pt, params)

        [k0, da, ~, a_ff, n_ff, tolabs, tolrel] = feval (@(x) x{:} , num2cell(params));        
        n = length(strt_pt);
        strt_pt(:,end+1) = strt_pt(:,1); %for ease while calculating Dr for the last starting point
        
        % The field at a point 'p' is given by:
        % phi(p) = phi_inc(p) + phi_scat(p)
        % phi_scat(p) = - oint[g1(r,p) grad(phi(r)).n - grad(g1(r,p).n phi(r)]dr (1)
        % where the closed loop integral is over the boundary, where we know the fields
        
        % 1. Computing the scattered field at the far field points
        scat_field_ff = zeros(1,n_ff);
        if eps_r == Inf %implies PEC
            for i = 1:n_ff
                for j = 1:n
                    Dr = strt_pt(:,j+1) - strt_pt(:,j);
                    % For PEC the second term in eqn (1) is zero. First term:
                    integral_g1 = da*integral(@(d)green2d(k0, ff_pt(:,i), strt_pt(:,j), Dr, d),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);
                    scat_field_ff(i) = scat_field_ff(i) - fields_bndry(j)*integral_g1;
                end
            end
        else
            for i = 1:n_ff
                for j = 1:n
                    Dr = strt_pt(:,j+1) - strt_pt(:,j);
                    t_hat = Dr/norm(Dr);
                    n_hat = [-t_hat(2), t_hat(1)];
                    integral_g1 = da*integral(@(d)green2d(k0, ff_pt(:,i), strt_pt(:,j), Dr, d),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);
                    integral_grad_g1_dot_n = da*integral(@(d)gradgreen2d_dot_n(k0, ff_pt(:,i), strt_pt(:,j), Dr, d, n_hat),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);
                    scat_field_ff(i) = scat_field_ff(i) - (fields_bndry(j)*integral_g1 - fields_bndry(j+n)*integral_grad_g1_dot_n);
                end
            end
        end

        RCS = 2*pi*a_ff*abs(scat_field_ff).^2; % (the magnitude of the incident field is just 1)
end
