E=200e3; L=1e3;  diam=40; a=40; P=1e3;
I_bar = pi*diam^4/64; I_beam = a^4/12; A_bar = pi*diam^4/4; A_beam=a^2;
Nnodes=9; Nelem=6; dofs=20;

elem =[1 2; 2 3; 4 5; 5 7;6 8; 3 9];

nodes= [0 0; 0.7071*L .7071*L ; 1.7071*L .7071*L ; 1.7071*L .7071*L;...
        1.353*L .353*L ;1.353*L .353*L ; L 0; 2*L 0; 2.7071*L .7071 *L];

BARRA = true;
esbarra=[BARRA 0 0 0 BARRA 0];
esbarra= logical(esbarra);
nodeDofs=[1 2 0; 3 4 5;6 7 8;6 7 9 ;10 11 12;...
          10 11 0; 13 14 15;16 17 0 ;18 19 20];
      
 %armo k_g
 k_g= zeros(dofs);
 for i=1:length(elem)
     if esbarra(i)
         ubic = [nodeDofs(elem(i,1),1:2),nodeDofs(elem(i,2),1:2)];
         L = norm( nodes(elem(i,2),:) - nodes(elem(i,1),:));
         k_loc = A_bar*E/L * [1 -1;...
                             -1 1];
                         
         dts = (nodes(elem(i,2),:)-nodes(elem(i,1),:))/L ;
         T=[dts 0 0; 0 0 dts];
         k_g(ubic,ubic) = T'*k_loc*T + k_g(ubic,ubic);
     else
         ubic = [nodeDofs(elem(i,1),:),nodeDofs(elem(i,2),:)];
         L = norm( nodes(elem(i,2),:) - nodes(elem(i,1),:));
         X = A_beam*E/L; Y_1=12*E*I_beam/L^3; Y_2=6*E*I_beam/L^2; Y_3=4*E*I_beam/L; Y_4 = 2*E*I_beam/L;  
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
 end
 
      
      
      
fzas = zeros(dofs,1); fzas(19) = -P;
bc= zeros(dofs,1); 
bc(1:2) = [1 1]; bc(13:14) = [1 1]; bc(16:17) = [1 1];


fzas_red= fzas(~bc);
k_red=k_g(~bc,~bc);

despl=zeros(dofs,1);
despl(~bc)=k_red^(-1)*fzas_red;

%pos-proceso
sigma_axial = E*[-1/L 1/L] *[despl(3) despl(6)]';
sigma_bending = E*[-6/L^2+12*L/2/L^3, -2/L+6*L/2/L^2, 6/L^2-12*L/2/L^3, -2/L+6*L/2/L^2]...
                 * [despl(4) despl(5) despl(7) despl(8)]';
sigma_b= @(y) sigma_axial+sigma_bending*y; 

%fplot(sigma_b,[-20,20]);
%view([90 -90]);

i=3;
ubic = [nodeDofs(elem(i,1),:),nodeDofs(elem(i,2),:)];
L = norm( nodes(elem(i,2),:) - nodes(elem(i,1),:));
X = A_beam*E/L; Y_1=12*E*I_beam/L^3; Y_2=6*E*I_beam/L^2; Y_3=4*E*I_beam/L; Y_4 = 2*E*I_beam/L;  
k_loc =  [X 0 0 -X 0 0;...
          0 Y_1 Y_2 0 -Y_1 Y_2;...
          0 Y_2 Y_3 0 -Y_2 Y_4;...
         -X 0 0 X 0 0;...
          0 -Y_1 -Y_2 0 Y_1 -Y_2;...
          0 Y_2 Y_4 0 -Y_2 Y_3];
                         
dts = (nodes(elem(i,2),:)-nodes(elem(i,1),:))/L ;
T=zeros(6); T(1,1:2)=dts; T(4,[4 5])=dts; T(2,1:2)= [-1*dts(2) dts(1)]; T(5,4:5)= [-1*dts(2) dts(1)];
T(3,3) =1;T(6,6)=1;

despl_3 = [despl([6 7 9]); despl(10:12)];
despl_local = T * despl_3; 

fzas_loc=k_loc*despl_local;
fzas_global = T'*fzas_loc;





