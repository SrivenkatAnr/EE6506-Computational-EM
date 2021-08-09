function [eps_arr] = get_multilayer_eps(seq, n, eps_n)
        %seq=1: alternate array, seq=2: fibonocci
        %n: number of terms in the sequence
        %assumes only 2 materials: material N and Air
        
        eps_arr = [];
        if seq == 1    %alternate arrangement
                for i = 1:n-1
                        eps_arr(end+1:end+2) = [eps_n 1];
                end
                eps_arr(end+1) = eps_n;         
        else    %fibonocci arrangement
                f1 = [1];
                f2 = [eps_n];
                if n == 1
                        eps_arr = f1;
                elseif n == 2
                        eps_arr = f2;
                else 
                        for i = 3:n
                                eps_arr = cat(2, f2, f1);
                                f1 = f2;
                                f2 = eps_arr;
                        end        
                end             
         end