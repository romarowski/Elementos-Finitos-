E=210e3; a=100; d=50;
Av=a^2; Iz=a^4/12; Ap=pi*d^2/4;

nodeDofs=[1:5 0 6:20];
nodeDofs=reshape(nodeDofs,3,7); nodeDofs=nodeDofs';

youngs=[E E 50*E];

elem=[2 3;3 6;4 7;1 3;3 4;5 6;6 7];
nodes=[0 0;62 -6;32 24; 72 60; 0 24; 32 48;72 84]; nodes=nodes*10;

k_g=zeros(20);
m=1;
for i=1:7
    if i<4
       ubic= [nodeDofs(elem(i,1),1:2) , nodeDofs(elem(i,2),1:2)];
       Le=norm(nodes(elem(i,2),:) - nodes(elem(i,1),:));
       X = youngs(m)*Ap/Le;
       k_loc = X*[1 -1;-1 1];
       dts = (nodes(elem(i,2),:) - nodes(elem(i,1),:))/Le;
       T=[dts 0 0;0 0 dts];
       k_g(ubic,ubic)= T'*k_loc*T + k_g(ubic,ubic);
       m=m+1;
    else
       ubic= [nodeDofs(elem(i,1),:) , nodeDofs(elem(i,2),:)];
       Le=norm(nodes(elem(i,2),:) - nodes(elem(i,1),:));
       X = E*Av/Le;
    Y1 = 12*E*Iz/Le^3;
    Y2 = 6*E*Iz/Le^2;
    Y3 = 4*E*Iz/Le^1;
    Y4 = 2*E*Iz/Le^1;
    
    Kdiag = diag([X Y1 Y3 X Y1 Y3]);
    Kp = [0     0       0       -X      0       0
          0     0       Y2      0       -Y1     Y2
          0     0       0       0       -Y2     Y4
          0     0       0       0       0       0
          0     0       0       0       0       -Y2
          0     0       0       0       0       0];
      k_loc = Kp + Kp' + Kdiag;

       
       dts = (nodes(elem(i,2),:) - nodes(elem(i,1),:))/Le;
       T=zeros(6);  T(1,1:2) =dts;T(4,4:5) =dts; 
       T(2,1:2) =[-1*dts(2) dts(1)]; T(5,4:5) =[-1*dts(2) dts(1)];
       T(3,3)=1; T(6,6)=1;
       k_g(ubic,ubic)= T'*k_loc*T + k_g(ubic,ubic);
    end
end

bc=zeros(20,1);
bc(1:3) =ones(3,1); bc(12:14) =ones(3,1); bc(4:5) =ones(2,1);

k_red=k_g(~bc,~bc);

despl=zeros(20,1);      

fzas=zeros(20,1);
fzas(19)=-4e3; fzas(10)=-4e3; fzas(18)=-9.375e3; fzas(19)=9.375e3;
fzas_red=fzas(~bc);

despl(~bc)=k_red^(-1)*fzas_red;

despl_e=despl(7)

i=4;
   ubic= [nodeDofs(elem(i,1),:) , nodeDofs(elem(i,2),:)];
       Le=norm(nodes(elem(i,2),:) - nodes(elem(i,1),:));
       X = E*Av/Le;
    Y1 = 12*E*Iz/Le^3;
    Y2 = 6*E*Iz/Le^2;
    Y3 = 4*E*Iz/Le^1;
    Y4 = 2*E*Iz/Le^1;
    
    Kdiag = diag([X Y1 Y3 X Y1 Y3]);
    Kp = [0     0       0       -X      0       0
          0     0       Y2      0       -Y1     Y2
          0     0       0       0       -Y2     Y4
          0     0       0       0       0       0
          0     0       0       0       0       -Y2
          0     0       0       0       0       0];
      k_loc = Kp + Kp' + Kdiag;

       
       dts = (nodes(elem(i,2),:) - nodes(elem(i,1),:))/Le;
       T=zeros(6);  T(1,1:2) =dts;T(4,4:5) =dts; 
       T(2,1:2) =[-1*dts(2) dts(1)]; T(5,4:5) =[-1*dts(2) dts(1)];
    T(3,3)=1; T(6,6)=1;
    despl_loc=T*despl(ubic);
    fzas_loc=k_loc*despl_loc;
    fzas_glob=zeros(2,6);
    fzas_glob(1,:)=T'*fzas_loc;
    
    
i=5;
   ubic= [nodeDofs(elem(i,1),:) , nodeDofs(elem(i,2),:)];
       Le=norm(nodes(elem(i,2),:) - nodes(elem(i,1),:));
       X = E*Av/Le;
    Y1 = 12*E*Iz/Le^3;
    Y2 = 6*E*Iz/Le^2;
    Y3 = 4*E*Iz/Le^1;
    Y4 = 2*E*Iz/Le^1;
    
    Kdiag = diag([X Y1 Y3 X Y1 Y3]);
    Kp = [0     0       0       -X      0       0
          0     0       Y2      0       -Y1     Y2
          0     0       0       0       -Y2     Y4
          0     0       0       0       0       0
          0     0       0       0       0       -Y2
          0     0       0       0       0       0];
      k_loc = Kp + Kp' + Kdiag;

       
       dts = (nodes(elem(i,2),:) - nodes(elem(i,1),:))/Le;
       T=zeros(6);  T(1,1:2) =dts;T(4,4:5) =dts; 
       T(2,1:2) =[-1*dts(2) dts(1)]; T(5,4:5) =[-1*dts(2) dts(1)];
    T(3,3)=1; T(6,6)=1;
    despl_loc=T*despl(ubic);
    fzas_loc=k_loc*despl_loc;
    
    fzas_glob(2,:)=T'*fzas_loc;
   
    i=2;
    
    ubic= [nodeDofs(elem(i,1),1:2) , nodeDofs(elem(i,2),1:2)];
       Le=norm(nodes(elem(i,2),:) - nodes(elem(i,1),:));
       X = E*Ap/Le;
       k_loc = X*[1 -1;-1 1];
       dts = (nodes(elem(i,2),:) - nodes(elem(i,1),:))/Le;
       T=[dts 0 0;0 0 dts];
    despl_loc=T*despl(ubic);
    fzas_loc=k_loc*despl_loc;
    
    fzas_bar=T'*fzas_loc;
       
    
    i=1;
    
    ubic= [nodeDofs(elem(i,1),1:2) , nodeDofs(elem(i,2),1:2)];
       Le=norm(nodes(elem(i,2),:) - nodes(elem(i,1),:));
       X = E*Ap/Le;
       k_loc = X*[1 -1;-1 1];
       dts = (nodes(elem(i,2),:) - nodes(elem(i,1),:))/Le;
       T=[dts 0 0;0 0 dts];
    despl_loc=T*despl(ubic);
    fzas_loc=k_loc*despl_loc;
    
    fzas_piston=T'*fzas_loc;
    fzas_piston=norm(fzas_piston(3:4))
       
    
    fzapist=fzas_glob(1,4:5)+fzas_glob(2,1:2)+fzas_bar(1:2)';
    
    i=7;
   ubic= [nodeDofs(elem(i,1),:) , nodeDofs(elem(i,2),:)];
       L=norm(nodes(elem(i,2),:) - nodes(elem(i,1),:));
       
       dts = (nodes(elem(i,2),:) - nodes(elem(i,1),:))/Le;
       T=zeros(6);  T(1,1:2) =dts;T(4,4:5) =dts; 
       T(2,1:2) =[-1*dts(2) dts(1)]; T(5,4:5) =[-1*dts(2) dts(1)];
    T(3,3)=1; T(6,6)=1;
    despl_loc=T*despl(ubic);
    
    B_bend=[-6/L^2+12*L/2/L^3 -4/L+6*L/2/L^2 6/L^2-12*L/2/L^3 -2/L+6*L/2/L^2];
    B_ax=[-1/L 1/L];
    sigma_ax= E*B_ax*despl_loc([1 4]);
    sigma_b= E*B_bend*despl_loc([2 3 5 6])*a/2;
    sigma_tot= sigma_ax-sigma_b;
    