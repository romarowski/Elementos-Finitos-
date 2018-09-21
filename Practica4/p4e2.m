E=210e3; A=1e3;

nodeDofs=1:15; nodeDofs = reshape(nodeDofs,3,5); nodeDofs=nodeDofs';
elem=[1 2;1 3; 1 4; 1 5]; dofs=15;
nodes=[4 4 3;0 4 0;0 4 6; 4 0 3; 8 -1 1];

k_g= zeros(dofs);
for i=1:length(elem)
ubic = [nodeDofs(elem(i,1),:),nodeDofs(elem(i,2),:)];
L = norm( nodes(elem(i,2),:) - nodes(elem(i,1),:));
k_loc = A*E/L * [1 -1;...
                -1 1];
                         
dts = (nodes(elem(i,2),:)-nodes(elem(i,1),:))/L ;
T=[dts 0 0 0; 0 0 0 dts];
k_g(ubic,ubic) = T'*k_loc*T + k_g(ubic,ubic);
end

bc=ones(dofs,1); bc(1:3)=[0 0 0];
fzas=zeros(dofs,1); %fzas(2)=-10e3;
c=[2 4:15]; x=[1 3];

fzas_red= fzas(~bc); k_red=k_g(~bc,~bc);

k_cc=k_g(c,c); k_xx = k_g(x,x); k_xc =k_g(x,c);

despl_c = zeros(length(c),1); despl_c(1) = -1.5171*.1;
despl_x=k_xx^(-1)*(-k_xc*despl_c);

despl=zeros(dofs,1);
despl(c)= despl_c; despl(x)=despl_x;

sigma_axial=zeros(4,1);
for i=1:length(elem)
ubic = [nodeDofs(elem(i,1),:),nodeDofs(elem(i,2),:)];
L = norm( nodes(elem(i,2),:) - nodes(elem(i,1),:));
dts = (nodes(elem(i,2),:)-nodes(elem(i,1),:))/L ;
T=[dts 0 0 0; 0 0 0 dts];
despl_local= T*despl(ubic);
sigma_axial(i) = E*[-1/L 1/L] *despl_local;
end

