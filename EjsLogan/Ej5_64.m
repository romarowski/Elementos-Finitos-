E=18e6; I=4.5e-4; A= 0.06;
Nelem=20; tita=90; R=2; Nnodes=Nelem+1;
dof_beam =3; dof_gral= dof_beam*Nnodes;

angulos =90:-180/Nelem:-90;
equis = R*cosd(angulos); yi = R*sind(angulos);
nodes=[equis' yi'];

R_new = 2.015915;
equis_new =R_new*cosd(angulos); yi_new = R_new*sind(angulos);
nodes_new = [equis_new' yi_new'];
despl_known=nodes_new-nodes;

elem = zeros(Nelem,2);
elem(:,1) = 1:Nelem; elem(:,2) = 2:Nnodes;

k_g = zeros(dof_gral);

for i=1:length(elem)
   ubic_gral=[3*elem(i,1)-2:3*elem(i,1) 3*elem(i,2)-2:3*elem(i,2)];
   dts= nodes(elem(i,2),:) -nodes(elem(i,1),:);
   largo = norm(dts);
   dts=dts/largo;
   phi=atan(dts(2)/dts(1));
   k_g(ubic_gral,ubic_gral) = k_g(ubic_gral,ubic_gral) +elemviga(E,I,A,largo,1,phi);
end


bc = false(3,Nnodes);
bc(3,:)=1; ; 
bc= reshape(bc,numel(bc),1);

k_uknown = k_g(~bc,~bc);
k_known = k_g(~bc,~bc);
k_xc = k_g(~bc,bc);

%despl = zeros(3,Nnodes);

%despl_uknown = k_uknown^(-1)*(-k_xc*despl_known);
despl_known = reshape(despl_known',numel(despl_known),1);
fzas=k_known*despl_known;
%fzas = zeros(3,Nnodes);





