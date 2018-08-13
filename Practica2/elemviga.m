function [k_e] = elemviga( E,I,A,L,N,phi)
% matriz en coordenadas genericas para elem viga doblado
syms X Y_1 Y_2 Y_3 Y_4;
dof = 3;
dof_g = (N+1)*dof;
X = A*E/(L/N);
Y_1 = 12*E*I/(L/N)^3;
Y_2 = 6*E*I/(L/N)^2;
Y_3 = 4*E*I/(L/N);
Y_4 = 2*E*I/(L/N);
k = [X 0 0 -X 0 0;
    0 Y_1 Y_2 0 -Y_1 Y_2;
    0 Y_2 Y_3 0 -Y_2 Y_4;
    -X 0 0 X 0 0;
    0 -Y_1 -Y_2 0 Y_1 -Y_2;
    0 Y_2 Y_4 0 -Y_2 Y_3];
T=zeros(dof*2);
subT = [cos(phi) sin(phi) 0;
        -sin(phi) cos(phi) 0;
        0 0 1];
T(1:3,1:3) = subT;
T(4:end, 4:end) = subT;
elem = zeros(N,2);
for i=1:N
    elem(i, :) = [i i+1];
end
k_e = zeros(dof_g);
for i=1:N
    k_2d = T'*k*T;
    ubicengral = 3*elem(i,1)-2:3*elem(i,2);
    k_e(ubicengral,ubicengral)  = k_2d + k_e(ubicengral,ubicengral);
end


end

