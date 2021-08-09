function [test_pt, strt_pt] = get_shape_coords(a, da, n_sides, theta)
        %generates a regular polygon in 2D, oriented along theta, centred at origin
        ang = 360/n_sides;    
        c = a/(2*sind(ang/2));     %length of line segment connecting centre to vertex
        vrtx_A = [c*cosd(theta); c*sind(theta)];    %vertex pointing along theta
        dir_AB = [-cosd(90 - (ang/2) - theta); sind(90 - (ang/2) - theta)];
        rotat_ang = [cosd(ang) -sind(ang); sind(ang) cosd(ang)];         
        n_e = round(a/da);    %total no. of discretisation points
 
        test_pt = zeros(2, n_sides*n_e);    %points of delta testing
        strt_pt = zeros(2, n_sides*n_e);     %starting points of pulse basis
        test_pt(:,1:n_e) = vrtx_A + [da/2:da:a].*dir_AB;
        strt_pt(:,1:n_e) = vrtx_A + [0:da:a-da].*dir_AB;        
        for i = 1+n_e : n_e:n_sides*n_e
            test_pt(:,i:(i+n_e-1)) = rotat_ang*test_pt(:,i-n_e:i-1);
            strt_pt(:,i:(i+n_e-1)) = rotat_ang*strt_pt(:,i-n_e:i-1);
        end
        
end 
        
