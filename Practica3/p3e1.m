E=210e3; Iz= 1.33e4; A= 400; load = 5e3;
dof_beam = 3; nod_total = 8; dof_gral = dof_beam*nod_total;

elem=[1 2; 2 3; 3 4; 4 5 ; 3 5 ;4 6; 6 7; 3 7;4 7; 7 8];
phi_elem = [90  0  0 -90 -45 90 180 90 135 180];
phi_elem = deg2rad(phi_elem);
largos= [1e3 0.5e3 1e3 1e3 1.41e3 1e3 1e3 1e3 1.41e3 .5e3];
youngs = [E/2 E E E E E E E E E ];

k_g = zeros(dof_gral);
for i=1:length(elem)
    ubic_gral = [3*elem(i,1)-2:3*elem(i,1) 3*elem(i,2)-2:3*elem(i,2)];
    k_g(ubic_gral,ubic_gral) = k_g(ubic_gral,ubic_gral) + elemviga(youngs(i),Iz,A,largos(i),1,phi_elem(i));
end

bc_sim = zeros(nod_total,dof_beam);
bc_sim(1,:) =[1 1 1]; bc_sim(5,:) = [1 1 1]; bc_sim(2,[1 3]) =[1 1]; bc_sim(8,[1 3])= [1 1];
bc_sim= reshape(bc_sim',[dof_gral 1]);

bc_asim = zeros(nod_total,dof_beam);
bc_asim(1,:) =[1 1 1]; bc_asim(5,:) = [1 1 1]; bc_asim(2,2) =1; bc_asim(8,2)= 1;
bc_asim= reshape(bc_asim',[dof_gral 1]);



fzas = zeros(nod_total,dof_beam);
fzas(8,2) = load; fzas(6,1) = load;
fzas = reshape(fzas',[dof_gral 1]);
fzas_red = fzas(~bc_sim);
fzas_red_asim = fzas(~bc_asim);

k_red_sim= k_g(~bc_sim,~bc_sim);
k_red_asim= k_g(~bc_asim,~bc_asim);

despl_sim= zeros(dof_gral,1); despl_asim= despl_sim;
despl_sim(~bc_sim) = k_red_sim^(-1)*fzas_red; 
despl_asim(~bc_asim) = k_red_asim^(-1)*fzas_red_asim;

despl_sim = reshape(despl_sim, [3 8]);
despl_asim = reshape(despl_asim, [3 8]);
