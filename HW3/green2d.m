function g = green2d(k, r_p, r0, Dr, d)
        % r_p = [r_p_x r_p_y] is the primed coordinate
        % r0 = [r0_x r0_y] is the starting point
        % Dr = [Dr_x Dr_y] is the total change in r0
        % d \in [0,1) gives r = r0 + d*Dr

        % Green's function, g(r,r_p) = -j/4*H_0^(2)(k|r' - r|)

        R_x = (r0(1) + d.*Dr(1)) - r_p(1);
        R_y = (r0(2) + d.*Dr(2)) - r_p(2);
        norm_R = sqrt( R_x.^2 + R_y.^2 );
        g = -1j/4*besselh(0, 2, k.*norm_R);
end