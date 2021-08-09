function grad_g_dot_n = gradgreen2d_dot_n(k, r_p, r0, Dr, d, n_hat)
        % n = [nx ny] is the direction of the normal
        % rest of the variables are same as in green2d

        % grad(g).n = jk/4 H_1^(2)(k|r'-r|)hat(r-r').n

        R_x = (r0(1) + d.*Dr(1)) - r_p(1);
        R_y = (r0(2) + d.*Dr(2)) - r_p(2);
        norm_R = sqrt( R_x.^2 + R_y.^2 );
        hat_R_dot_n = (n_hat(1).*R_x + n_hat(2).*R_y)./norm_R;

        grad_g_dot_n = 1j*k/4*besselh(1, 2, k.*norm_R).*hat_R_dot_n;
end