syms x y z w mu;    %define symbols for x,y,z, permiability

E = [x * sin(z), y^2, z^3 * x];    %give input electric field expression here
%mu = 4*pi*(10^(-7));    %give permeability value here

H = symbolic_H_finder(E, mu)

function H = symbolic_H_finder(E, mu)
    syms x y z w;
    Ex = E(1); Ey = E(2); Ez = E(3);    %define x,y,z components of E
    
    %curl(E) = -jwu*H,    H = (j/wu) * curl(E)
    Hx = ( diff(Ez,y) - diff(Ey,z) ) * 1i / (w * mu) ;    %assign corresponding curl elements to components of H
    Hy = ( diff(Ex,z) - diff(Ez,x) ) * 1i / (w * mu) ;
    Hz = ( diff(Ey,x) - diff(Ex,y) ) * 1i / (w * mu) ;
    
    H = [Hx, Hy, Hz];
end    