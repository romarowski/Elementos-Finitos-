Largo=1e3; E=210e3; a=15; k=100; q=1;
A_beam=a^2; I_beam=a^4/12;


nodeDofs=1:18; nodeDofs=reshape(nodeDofs,3,6); nodeDofs=nodeDofs';
nodeDofs=[nodeDofs;19 20 0;21 22 0];
dofs=22;

elem=[1 2;2 3;3 4;4 5;3 6;1 7;5 8];
isBarra=[0 0 0 0 0 1 1]; isBarra=logical(isBarra);
dts_spring=[0 1;0 -1];

%r=0.1;
for r=.1:.01:.9 %abs(despl(2)) > 1e-10 || abs(despl(14)) > 1e-10 

nodes=[0 0; r*Largo/2 0;Largo/2 0;Largo-r*Largo/2 0;Largo 0;Largo/2 Largo/2];

k_g=zeros(dofs); j=1; m=1;
fzas=zeros(dofs,1);
for i=1:length(elem)
     if isBarra(i)
         ubic = [nodeDofs(elem(i,1),1:2),nodeDofs(elem(i,2),1:2)];
         k_loc = k * [1 -1;...
                      -1 1];
                         
         dts = dts_spring(j,:) ;
         T=[dts 0 0; 0 0 dts];
         k_g(ubic,ubic) = T'*k_loc*T + k_g(ubic,ubic);
         j=j+1;
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
         %cargas
         if i<5
         fzas(m:m+5)=fzas(m:m+5)+ [0,-q*L/2,-q*L^2/12,0,-q*L/2,q*L^2/12]';
         m=m+3;
         end
     end
end
 
bc=zeros(dofs,1); bc(5)=1; bc(7:8)=[1 1]; bc(11) =1; bc(16:18) = [1 1 1];
bc(19:20)=[1 1]; bc(21:22)=[1 1];

k_red=k_g(~bc,~bc);

fzas_red=fzas(~bc);

despl=zeros(dofs,1);

despl(~bc)=k_red^(-1)*fzas_red;
% %r=r+0.01;


if abs(despl(2))<1e-2 && abs(despl(14))< 1e-2
    r
end
end

