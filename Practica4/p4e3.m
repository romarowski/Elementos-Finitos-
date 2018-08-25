E =210e3; I_beam=2e8; A=2e4;
dofs=12;

nodeDofs= 1:12; nodeDofs=reshape(nodeDofs',[3 4]); nodeDofs= nodeDofs';
nodes=[0 0; 0 6e3; 6e3  6e3; 6e3 0];
elem=[1 2;2 3;3 4];

 k_g= zeros(dofs);
 for i=1:length(elem)
         ubic = [nodeDofs(elem(i,1),:),nodeDofs(elem(i,2),:)];
         L = norm( nodes(elem(i,2),:) - nodes(elem(i,1),:));
         X = A*E/L; Y_1=12*E*I_beam/L^3; Y_2=6*E*I_beam/L^2; Y_3=4*E*I_beam/L; Y_4 = 2*E*I_beam/L;  
         k_loc =  [X 0 0 -X 0 0;...
                   0 Y_1 Y_2 0 -Y_1 Y_2;...
                   0 Y_2 Y_3 0 -Y_2 Y_4;...
                   -X 0 0 X 0 0;...
                   0 -Y_1 -Y_2 0 Y_1 -Y_2;...
                   0 Y_2 Y_4 0 -Y_2 Y_3];
                         
         dts = (nodes(elem(i,2),:)-nodes(elem(i,1),:))/L ;
         T=zeros(6); T(1,1:2)=dts; T(4,[4 5])=dts; T(2,1:2)= [-1*dts(2) dts(1)]; T(5,4:5)= [-1*dts(2) dts(1)];
         T(3,3) =1;T(6,6)=1;
         k_g(ubic,ubic) = T'*k_loc*T + k_g(ubic,ubic); 
 end
 
 bc=zeros(dofs,1); bc(1:3)= ones(3,1); bc(10:12)= ones(3,1); 
 
 fzas=zeros(dofs,1); fzas(4)=15e3; fzas(6)=10e6;
 fzas_red=fzas(~bc);
 
k_red=k_g(~bc,~bc);

despl= zeros(dofs,1);

despl(~bc) = k_red^(-1)*fzas_red;

despl_loc=[0.0071 4.2871 -0.0005 0 0 0];

sigma_axial = E*[-1/6e3 1/6e3]*[despl_loc(1) 0]';
sigma_bending = E * [-6/L^2+12/L^2 -4/L+6/L 6/L^2-12/L^2 -2/L+6/L]*despl_loc(2:5)';

sigma_tot =sigma_bending+sigma_axial;