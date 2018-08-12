syms E A L13 L24 I phi;
E = 210e3;
A = .5e-2*1000^2;
L1= 8e3;
L2 = 10e3;
I = .5e-4*1000^4;
P = 300;
dof = 3;

Nvert = 100; %debe ser par! 
Nhor = 1;
N_total = Nvert + 2*Nhor;

k_evert = elemviga(E, I, A, L1/Nvert,Nvert,deg2rad(90));
k_ehor = elemviga(E, I, A, L2,Nhor,0);
dof_g = (N_total +1)*dof;

nodes = [zeros(Nvert+1,1) [0:L1/Nvert:L1]'];

elem = zeros(3,2);
elem(1,:) = [1 Nvert+1];
node_a_4000 = ceil(length(nodes)/2); 
elem(2,:) = [node_a_4000 Nvert+3];
ult_node_vert = Nvert+1;
elem(3,:) = [Nvert+1 Nvert+2];

%k_total
k_total = zeros(dof_g);

ubicengral = 3*elem(1,1)-2:3*elem(1,2); %ensamble elementos verticales
k_total(ubicengral,ubicengral)  = k_evert + k_total(ubicengral,ubicengral);

ubicengral_1 = 3*elem(2,1)-2:3*elem(2,1);        %ensamble elemento horizontal abajo
ubicengral_2 = 3*elem(2,2)-2:3*elem(2,2); 
k_total(ubicengral_1,ubicengral_1)  = k_ehor(1:3,1:3) + k_total(ubicengral_1,ubicengral_1);      
k_total(ubicengral_1,ubicengral_2)  = k_ehor(1:3,4:6) + k_total(ubicengral_1,ubicengral_2);      
k_total(ubicengral_2,ubicengral_1)  = k_ehor(4:6,1:3) + k_total(ubicengral_2,ubicengral_1);
k_total(ubicengral_2,ubicengral_2)  = k_ehor(4:6,4:6) + k_total(ubicengral_2,ubicengral_2);

ubicengral = 3*elem(3,1)-2:3*elem(3,2); %ensamble elemento horizontal arriba
k_total(ubicengral,ubicengral)  = k_ehor + k_total(ubicengral,ubicengral);


%vectorfzas
distload = @(y) P * (L1 - y)/ (L1);
dof_nodal = 3;
dof_elemental = 2*dof_nodal;
dof_vert = 3*(Nvert+1);
fzas = zeros(dof_g,1);
j=1;
for i=1:3:(dof_vert -5)
    L_e = (L1/Nvert);
    fzahor_sq = distload(nodes(j+1,2)) * L_e  /2  ;
    momento_sq = distload(nodes(j+1,2)) * (L_e)^2/12;
    fza_elem_sq =  [fzahor_sq 0 -momento_sq fzahor_sq 0 momento_sq]';
    
    q_triang = distload(nodes(j,2))-distload(nodes(j+1,2));
    fzahor_triang = q_triang*L_e;
    mom_triang = q_triang * L_e^2;
    fza_elem_tri = [7/20 * fzahor_triang, 0, -1/20 * mom_triang, 3/20 * fzahor_triang 0 mom_triang/30]';
    
    fzas(i:i+5) = fzas(i:i+5) + fza_elem_sq + fza_elem_tri;    
    j=j+1;
end


bc = zeros(dof_g,1);
bc(1:3) = ones(3,1);
bc((end-5):end) = ones(6,1);
k_red=k_total(~bc,~bc);
fzas_red= fzas(~bc);
despl=k_red^(-1)*fzas_red;
