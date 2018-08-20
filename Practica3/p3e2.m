E=210e3; A=1e4; load= 1250;
dof_bar =3; n_nodes =4; dof_gral =dof_bar*n_nodes;

elem=[1 2; 1 3; 2 4 ;1 4];
youngs=[E E E/2 E];
nodes= [10e3 0 0; 0 0 17.32e3; 0 0 0; 0 20e3 17.32e3];

k_g = zeros(dof_gral);

for i=1:length(elem)
   ubic_gral=[3*elem(i,1)-2:3*elem(i,1) 3*elem(i,2)-2:3*elem(i,2)];
   dts= nodes(elem(i,2),:) -nodes(elem(i,1),:);
   largo = norm(dts);
   dts=dts/largo;
   T = [dts 0 0 0
       0 0 0 dts];
   k_loc = A*youngs(i)/largo * [ 1 -1;-1 1];
   k_e= T'*k_loc*T;
   k_g(ubic_gral,ubic_gral) = k_g(ubic_gral,ubic_gral) +k_e;
end

bc = zeros(3,n_nodes);
bc(:,1) = [1 1 1]; bc(:,2) = [1 1 1]; bc(:,3) = [1 1 1]; bc(1,4) =1;
bc=reshape(bc,numel(bc),1);

fzas= zeros(3,n_nodes);
fzas(3,4) = -load;
fzas=reshape(fzas,numel(fzas),1);
fzas_red=fzas(~bc);

k_red= k_g(~bc,~bc);

despl= zeros(dof_gral,1);
despl(~bc)=k_red^(-1)*fzas_red;
despl_lindo=reshape(despl,[3 n_nodes]);

sigma=zeros(length(elem),1);
for i=1:length(elem)
   ubic_gral=[3*elem(i,1)-2:3*elem(i,1) 3*elem(i,2)-2:3*elem(i,2)];
   dts= nodes(elem(i,2),:) -nodes(elem(i,1),:);
   largo = norm(dts);
   dts=dts/largo;
   T = [dts 0 0 0
       0 0 0 dts];
   despl_loc= T*despl(ubic_gral);
   B=[-1/largo 1/largo];
   sigma(i) = E * B *despl_loc;
   
end


