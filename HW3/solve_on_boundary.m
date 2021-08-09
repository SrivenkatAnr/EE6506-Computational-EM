function x = solve_on_boundary(eps_r, test_pt, strt_pt, params, a_sc)
        %finds \phi and \grad \phi. \hat{n} at the boundary using Extinction theorem
        [k0, da, theta_i, ~, ~, tolabs, tolrel] = feval (@(x) x{:} , num2cell(params));         
        n = length(test_pt);
        strt_pt(:,end+1) = strt_pt(:,1); %for ease while calculating Dr for the last starting point
        
        %framing SI equations in Ax=b form
        if eps_r == Inf %implies PEC
            A = zeros(n, n);
            b = zeros(n,1);
            for i = 1:n
                for j = 1:n
                    Dr = strt_pt(:,j+1) - strt_pt(:,j);
                    A(i, j) = integral(@(d)green2d(k0, test_pt(:,i), strt_pt(:,j), Dr, d),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);
                end
                b(i) = inc_field(test_pt(:,i), theta_i, k0);
            end
            A = da*A;    %scaling by the length of the segment since the integral iterated over d from 0 to 1
            x = A\b;       %solve by gaussian elimination
        else
            A = zeros(2*n, 2*n);
            b = zeros(2*n,1);
            for i = 1:n
                for j = 1:n
                    Dr = strt_pt(:,j+1) - strt_pt(:,j);
                    t_hat = Dr/norm(Dr);
                    n_hat = [-t_hat(2), t_hat(1)];
                    if i==j    %considering cases of singularity separately
                         A(i, j)  = integral(@(d)green2d(k0, test_pt(:,i), strt_pt(:,j), Dr, d),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);           %The g1 term (integrable singularity)
                         A(i+n, j)= integral(@(d)green2d(sqrt(eps_r)*k0, test_pt(:,i), strt_pt(:,j), Dr, d),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);         %The g2 term (integrable singularity)
                         A(i,j+n) = 0.5 -1*integral(@(d)gradgreen2d_dot_n(k0, test_pt(:,i), strt_pt(:,j), Dr, d, n_hat),0.0,0.5-a_sc,'AbsTol',tolabs,'RelTol',tolrel)...
                             -1*integral(@(d)gradgreen2d_dot_n(k0, test_pt(:,i), strt_pt(:,j), Dr, d, n_hat),0.5+a_sc,1,'AbsTol',tolabs,'RelTol',tolrel);          %The grad(g1).n term 
                         A(i+n,j+n) = -0.5 -1*integral(@(d)gradgreen2d_dot_n(sqrt(eps_r)*k0, test_pt(:,i), strt_pt(:,j), Dr, d, n_hat),0.0,0.5-a_sc,'AbsTol',tolabs,'RelTol',tolrel)...
                             -1*integral(@(d)gradgreen2d_dot_n(sqrt(eps_r)*k0, test_pt(:,i), strt_pt(:,j), Dr, d, n_hat),0.5+a_sc,1,'AbsTol',tolabs,'RelTol',tolrel);     %The grad(g2).n term
                    else
                         A(i, j)     = integral(@(d)green2d(k0, test_pt(:,i), strt_pt(:,j), Dr, d),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);            %The g1 term
                         A(i+n, j)   = integral(@(d)green2d(sqrt(eps_r)*k0, test_pt(:,i), strt_pt(:,j), Dr, d),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);         %The g2 term 
                         A(i, j+n)   = -1*integral(@(d)gradgreen2d_dot_n(k0, test_pt(:,i), strt_pt(:,j), Dr, d, n_hat),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);            %The grad(g1).n term
                         A(i+n, j+n) = -1*integral(@(d)gradgreen2d_dot_n(sqrt(eps_r)*k0, test_pt(:,i), strt_pt(:,j), Dr, d, n_hat),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);      %The grad(g2).n term
                    end
                end
                b(i) = inc_field(test_pt(:,i), theta_i, k0);
            end
            A = da*A;    %scaling by the length of the segment since the integral iterated over d from 0 to 1
            x = A\b;       %solve by gaussian elimination
        end
end