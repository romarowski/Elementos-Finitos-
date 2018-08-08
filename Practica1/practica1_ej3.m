syms E A L k phi;
E = 210e3;
A = 5e4;
L = 5e3;
k = 4000;
P = -100e3/2;
k_spring = 4000;
phi = deg2rad(60);
N=1;
dof = 2;
dof_g = (N+1)*dof;
nodes = [0:-L*cos(phi)/N:-L*cos(phi); 0:L*sin(phi)/N:L*sin(phi)];
nodes = nodes';
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
bc= zeros(1,dof_g);
bc(end) = 1;
bc(end-1) =1;
bc(1)=1;
fzas = zeros(dof_g,1);
fzas(2) = P;
k_g(2,2) = k_g(2,2) + k_spring/2;
k_red = k_g(~bc,~bc);
despl = k_red^(-1)*fzas(2:end-2);
despl= [0 despl 0 0];


