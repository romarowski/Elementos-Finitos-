E=30457924.9; Apuente = 52.3; Asosten = 52.3; Abarra = 9.62; Lsosten= 600; Lpuente = 2400; Isosten = 6990; 
Ipuente = 6990; Lsuperpalo = 799.9; q_distr = 125; Ncables = 1; % la cantidad de cables es simetrica
dof_bar=2; dof_beam=3; Nelem_puente= 2*(Ncables+1); Nelem_superpalo = Ncables; Nsosten=1;
Nelem_total = Nelem_puente + Nelem_superpalo   + 2*Ncables + Nsosten;

elem=zeros(Nelem_total,2);
for i=1:Nelem_puente
    elem(i,:) = [i i+1];
end
node_mitad_puente = Nelem_puente/2 +1 ;
elem(Nelem_puente+1,:) = [node_mitad_puente Nelem_puente+3]; % +3 por el nodo al final y el nodo de abajo
for j=Nelem_puente+3:Nelem_puente+1+Nelem_superpalo
    elem(j,:) = [j j+1];
end
elem_barras_puente = [2:node_mitad_puente-1 node_mitad_puente+1:Nelem_puente];
elem_barras_superpalo =  Nelem_puente+2+Nelem_superpalo:-1:Nelem_puente+3;
elem_barras = [elem_barras_puente ;[elem_barras_superpalo flip(elem_barras_superpalo)]];

elem(Nelem_puente+1+Nelem_superpalo:end-1,:) = elem_barras';
                



