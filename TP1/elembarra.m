function [k_g] = elembarra(N,E, A, L, phi)
%  matriz de rigidez para elemento barra en el sist global
dof = 2;
dof_g = (N+1)*dof;

elem = zeros(N,2);
for i=1:N
    elem(i, :) = [i i+1];
end
transf = [cos(phi) sin(phi) 0 0 ;0 0 cos(phi) sin(phi)];
k_g = zeros(dof_g);
for i=1:N
    k_e = transf'* A * E / (L/N) * [1 -1 ; -1 1] * transf;
    ubicengral = 2*elem(i,1)-1:2*elem(i,2);
    k_g(ubicengral,ubicengral)  = k_e + k_g(ubicengral,ubicengral);
end

end

