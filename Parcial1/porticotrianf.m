E=210e3; Iz=.5e-4*1000^4; A=.5e-2*1000^2;

nodeDofs=1:18; nodeDofs=reshape(nodeDofs,3,6); nodeDofs=nodeDofs';
elem=[1 2;2 3;2 4;4 5;5 6];
nodes=[0 0;0 4;10 4;0 6; 0 8; 10 8]; nodes=nodes*1000;

k_g=zeros(18);

for i=1:5
    ubic=[nodeDofs(elem(i,1),:) nodeDofs(elem(i,2),:)];    
    Le= norm( nodes(elem(i,2),:) -nodes(elem(i,1),:));
    X = E*A/Le;
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
    Kel = Kp + Kp' + Kdiag;

    dts =  ( nodes(elem(i,2), : ) - nodes(elem(i,1),:) ) /Le;
    T=zeros(6); T(1,1:2) =dts; T(4,4:5) =dts; 
    T(2,1:2) =[-1*dts(2) dts(1)]; T(5,4:5) =[-1*dts(2) dts(1)];
    T(3,3)=1; T(6,6)=1;
    
    k_g(ubic,ubic)= T'*Kel*T +k_g(ubic,ubic);
end

bc=zeros(18,1);
bc(1:3)=ones(3,1); bc(7:9)=ones(3,1); bc(16:18)=ones(3,1);

k_red=k_g(~bc,~bc);

fzas=zeros(18,1);

fzas(4)=9e4+3e5+5.25e4+75e3; fzas(6) =8e7+2e8-15e6-25e6;
fzas(10)=22.5e3+75e3+5.25e4; fzas(12) =400e3+25e6-15e6;
fzas(13)=22.5e3; fzas(15)=400e3;

fzas_red=fzas(~bc);

despl=zeros(18,1);
despl(~bc)=k_red^-1*fzas_red;


i=3;
ubic=[nodeDofs(elem(i,1),:) nodeDofs(elem(i,2),:)];    
Le= norm( nodes(elem(i,2),:) -nodes(elem(i,1),:));
    X = E*A/Le;
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
    Kel = Kp + Kp' + Kdiag;

    dts =  ( nodes(elem(i,2), : ) - nodes(elem(i,1),:) ) /Le;
    T=zeros(6); T(1,1:2) =dts; T(4,4:5) =dts; 
    T(2,1:2) =[-1*dts(2) dts(1)]; T(5,4:5) =[-1*dts(2) dts(1)];
    T(3,3)=1; T(6,6)=1;
    
    despl_loc=T*despl(ubic);
    
    fzas_loc=Kel*despl_loc;
    fzas_en4=T'*fzas_loc;
    
