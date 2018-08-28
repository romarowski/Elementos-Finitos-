a=40; b=120; E=210e3; alfa=1e-6; dT=100;
A=a*b; Iz=a*b^3/12;

nodeDofs= [1 2 0;3 4 5;6 7 8;9 10 11; 12 13 14];
elem=[1 2;2 3;3 4;4 5]; 
nodes=[0 0;0 280;112.58 345;138.56 360;199.19 395];

isBarra=logical([1 0 0 0]);

dofs=14;

k_g=zeros(dofs);
for i=1:length(elem)
    if isBarra(i)
         ubic = [nodeDofs(elem(i,1),1:2),nodeDofs(elem(i,2),1:2)];
         L = norm( nodes(elem(i,2),:) - nodes(elem(i,1),:));
         k_loc = A*E/L * [1 -1;...
                             -1 1];
                         
         dts = (nodes(elem(i,2),:)-nodes(elem(i,1),:))/L ;
         T=[dts 0 0; 0 0 dts];
         k_g(ubic,ubic) = T'*k_loc*T + k_g(ubic,ubic);
    else    
        ubic = [nodeDofs(elem(i,1),:),nodeDofs(elem(i,2),:)];
         L = norm( nodes(elem(i,2),:) - nodes(elem(i,1),:));
         X = A*E/L; Y_1=12*E*Iz/L^3; Y_2=6*E*Iz/L^2; Y_3=4*E*Iz/L; Y_4 = 2*E*Iz/L;  
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

bc=zeros(dofs,1);
bc(1:2)=[1 1]; bc(12:14)= [1 1 1];

k_red=k_g(~bc,~bc);