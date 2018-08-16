Eacero=3e7; Econcreto= 64e5; Aviga = 62.2e2; Asosten = pi*40^2/4; Abarra = pi*5^2/4;
Lsosten= 300; Lpuente = 2400; Isosten = pi*40^4/64; Iviga=10300e4;
q_distr = -125; Ncables = 2; % la cantidad de cables es simetrica
dof_bar=2; dof_beam=3; Nelem_puente= 2*(Ncables+1); Nelem_superpalo = Ncables; Nsosten=1;
Nelem_total = Nelem_puente + Nelem_superpalo   + 2*Ncables + Nsosten;
nodos_totales = Nelem_superpalo + Nsosten + Nelem_puente +1; dof_total = nodos_totales*dof_beam;
nodes_puente = Nelem_puente +1; dof_puente = dof_beam*nodes_puente;
nodo_de_abajo= nodes_puente +1; 
seccion_viga_alto=5; %medida desde el centro

elem=zeros(Nelem_total,2);
for i=1:Nelem_puente
    elem(i,:) = [i i+1]; %cargo los del puente
end
node_mitad_puente = Nelem_puente/2 +1 ;
elem(Nelem_puente+1,:) = [node_mitad_puente Nelem_puente+3]; % +3 por el nodo al final y el nodo de abajo
for j=Nelem_puente+3:Nelem_puente+1+Nelem_superpalo
    elem(j-1,:) = [j j+1];
end
elem(end,:) = [node_mitad_puente i+2];
elem_barras_puente = [2:node_mitad_puente-1 node_mitad_puente+1:Nelem_puente];
elem_barras_superpalo =  Nelem_puente+2+Nelem_superpalo:-1:Nelem_puente+3;
elem_barras = [elem_barras_puente ;[elem_barras_superpalo flip(elem_barras_superpalo)]];
elem(Nelem_puente+1+Nelem_superpalo:end-1,:) = elem_barras';

%matriz de posicion de nodos
nodes = zeros (nodos_totales,2);
d_entre_cables_puente =  Lpuente / 2 / (Ncables+1);
d_entre_cables_superpalo = 0.1* d_entre_cables_puente;
nodos_horizontales = Lpuente/2:-d_entre_cables_puente:0;
nodos_verticales = nodos_horizontales(2):-d_entre_cables_superpalo:nodos_horizontales(2)-d_entre_cables_superpalo*(Ncables-1);

nodes(1:node_mitad_puente,1) = (-1)*nodos_horizontales;
nodes(node_mitad_puente:2*node_mitad_puente-1,1) = flip(nodos_horizontales);
nodes(2*node_mitad_puente+1:2*node_mitad_puente+1+ length(nodos_verticales)-1,2) = flip(nodos_verticales); 
nodes(2*node_mitad_puente,:) = [0 -300];

%sumo las vigas del puente
k_g = zeros(dof_total);
for i=1:Nelem_puente
    ubic_global = [3*elem(i,1)-2:3*elem(i,1)  3*elem(i,2)-2:3*elem(i,2) ];
    L = d_entre_cables_puente;
    k_g(ubic_global,ubic_global) = elemviga(Eacero,Iviga,Aviga,L,1,deg2rad(0)) + k_g(ubic_global,ubic_global); 
end
%sumo las vigas del superpalo
for j=Nelem_puente+1:Nelem_puente+Nelem_superpalo
    ubic_global = [3*elem(j,1)-2:3*elem(j,1)  3*elem(j,2)-2:3*elem(j,2) ];
    L = abs(nodes(elem(j,1),2)-nodes(elem(j,2),2));
    k_g(ubic_global,ubic_global) = elemviga(Econcreto,Isosten,Asosten,L,1,deg2rad(90)) + k_g(ubic_global,ubic_global); 
end
%sumo las barras
for j=Nelem_puente+Nelem_superpalo+1:Nelem_puente+Nelem_superpalo+2*Ncables
ubic_global = [3*elem(j,1)-2:3*elem(j,1)-1  3*elem(j,2)-2:3*elem(j,2)-1 ];
L = sqrt(nodes(elem(j,1),1)^2 + nodes(elem(j,2),2)^2);
if elem(j,1)<node_mitad_puente
    phi = asin(abs(nodes(elem(j,2),2)) / L);
else
   phi = pi/2+(pi/2-asin(nodes(elem(j,2),2) / L)) ;   
end
k_g(ubic_global,ubic_global) = elembarra(1,Eacero,Abarra,L,phi) + k_g(ubic_global,ubic_global);
end
%sumo la viga de abajo
ubic_global = [3*elem(end,1)-2:3*elem(end,1)  3*elem(end,2)-2:3*elem(end,2) ];
k_g(ubic_global,ubic_global) = elemviga(Econcreto,Isosten,Asosten,Lsosten,1,deg2rad(-90)) + k_g(ubic_global,ubic_global);

bc = zeros(dof_total,1);
bc(1:2) = [1 1] ; bc(dof_puente-1)=1; bc(nodo_de_abajo*dof_beam-2:nodo_de_abajo*dof_beam) = [1 1 1];

fzas = zeros(dof_total,1);
L_e_puente = d_entre_cables_puente; fza_vert = q_distr*L_e_puente/2; momento = q_distr*L_e_puente^2/12;
fza_elem_puente = [0 fza_vert momento 0 fza_vert -momento];
for i=1:3:3*nodes_puente-5
    fzas(i:i+5) = fzas(i:i+5) + fza_elem_puente';
end

fzas_red = fzas(~bc);
k_red = k_g(~bc,~bc);

despl = k_red^(-1)*fzas_red;

%post-proceso
despl = [[0 0] despl']; 
despl = [despl(1:nodes_puente*dof_beam-2) 0 despl(dof_beam*nodes_puente-1:end)];
despl = [despl(1:(nodes_puente+1)*dof_beam-3) [0 0 0] despl(dof_beam*(nodes_puente+1)-2:end)];

%tensiones normales en las vigas del puente
B_axial_puente = [-1/L_e_puente 1/L_e_puente];
despl_en_x = despl(1:3:dof_beam*nodes_puente);
sigma_axial_viga_puente =zeros(Nelem_puente,1);
for i=1:Nelem_puente
sigma_axial_viga_puente(i) = Eacero*B_axial_puente * despl_en_x([i i+1])';
end
%tensiones por flexion en las vigas del puente
despl_en_y_tita_viga_puente = zeros(nodes_puente*2,1);
n=1;
for i=2:3:dof_beam*nodes_puente
despl_en_y_tita_viga_puente([n n+1]) = despl([i i+1]);
n=n+2;
end
N_analisis_por_viga =5; sigma_bending_puente = zeros(Nelem_puente,N_analisis_por_viga+1);
n=1;
for i=1:Nelem_puente
    j=1;
    for x=0:L_e_puente/N_analisis_por_viga:L_e_puente
        y=seccion_viga_alto;
        B_flexion_puente = [-6/L_e_puente^2+12*x/L_e_puente^3; -4/L_e_puente+6*x/L_e_puente^2;
                            6/L_e_puente^2-12*x/L^3; -2/L+6*x/L^2];
        sigma_bending_puente(i,j) = y*Eacero*B_flexion_puente'*despl_en_y_tita_viga_puente(n:n+3);                       
        j=j+1;
    end
    n=n+2;
end
sigma_total_puente = zeros(size(sigma_bending_puente));
for i=1:Nelem_puente
    sigma_total_puente(i,:) = abs(sigma_bending_puente(i,:)) + abs(sigma_axial_viga_puente(i));
end

%tensiones en los cables

for j=Nelem_puente+Nelem_superpalo+1:Nelem_puente+Nelem_superpalo+2*Ncables
    L = sqrt(nodes(elem(j,1),1)^2 + nodes(elem(j,2),2)^2);
    if elem(j,1)<node_mitad_puente
        beta = asin(abs(nodes(elem(j,2),2)) / L);
    else
        beta = pi/2+(pi/2-asin(nodes(elem(j,2),2) / L)) ;   
    end
    transf = [cos(beta) sin(beta) 0 0 ;0 0 cos(beta) sin(beta)];
   
end
 
    

