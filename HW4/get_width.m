function [wid_arr] = get_width(eps_arr,wid, ratio)
        %generates array of widths of each layer
        
         wid_a = wid;
         wid_n = wid * ratio;
         
         wid_arr = zeros(size(eps_arr));
         for i = 1 : size(eps_arr')
                 if eps_arr(i) == 1 
                         wid_arr(i) = wid_a;    %set width of air gap
                 else
                         wid_arr(i) = wid_n;    %set width of material layer
                 end
         end        
         