syms c I E q P L1 L2;
I = 700e-6*1000^4;
E = 210e3;
q = -10;
P = -100e3;
N = 3;
L1 = 12e3;
L2 = 18e3;
dof = 2;
dof_g = (N+2)*dof;
nodes = 0:L1/N:L1; 
nodes = [nodes L2];
elem = zeros(N+1,2);
for i=1:N+1
        elem(i,:) = [i i+1];
end
k_g = zeros(dof_g);
for i=1:N+1
    if i~=N+1
        k_e =   E*I / (L1/N)^3 * [12 6*(L1/N) -12  6*(L1/N); 6*(L1/N) 4*(L1/N)^2 -6*(L1/N) 2*(L1/N)^2;-12 -6*(L1/N) 12  -6*(L1/N);6*(L1/N) 2*(L1/N)^2 -6*(L1/N) 4*(L1/N)^2 ];
    else
        k_e = E*2*I / (L1/N)^3 * [12 6*(L1/N) -12  6*(L1/N); 6*(L1/N) 4*(L1/N)^2 -6*(L1/N) 2*(L1/N)^2;-12 -6*(L1/N) 12  -6*(L1/N);6*(L1/N) 2*(L1/N)^2 -6*(L1/N) 4*(L1/N)^2 ];
    end
    ubicengral = 2*elem(i,1)-1:2*elem(i,2);
    k_g(ubicengral,ubicengral)  = k_e + k_g(ubicengral,ubicengral);
end
bc = zeros(1,dof_g);
bc(1) = 1;
bc(2) = 1;
bc(2*N+1) = 1;
fzas = zeros(dof_g,1);
fzas(end-1) = P;
fzavert = q * L1 / N /2;
momento = q * (L1/N)^2/12;
fzaelem = [fzavert momento fzavert -momento]';
for i=3:2:dof_g-5
    fzas(i:i+3) =  fzas(i:i+3) +fzaelem;
end
k_red = k_g(~bc,~bc);
despl = k_red^(-1)*fzas(~bc);

