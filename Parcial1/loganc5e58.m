E=30e6; G=10e3; A=100; Iz=1e3; Iy=200; Ix=100;

nodeDofs=[1:6;7:12;13:18;19:24];
elem=[1 2;2 3;3 4];
nodes=[0 0 0; 120 0 0;120 0 -120;120 -240 -120];
p_aux=[60 0 -10; 120 10 -60; 120 -120 -110];

k_g=zeros(24);
for i=1:3
    ubic=[nodeDofs(elem(i,1),:) nodeDofs(elem(i,2),:)];
    dir=nodes(elem(i,2),:)- nodes(elem(i,1),:);
    L=norm(nodes(elem(i,2),:)- nodes(elem(i,1),:));
    X=A*E/L;Y1=12*E*Iz/L^3;Z1=12*E*Iy/L^3;Y2=6*E*Iz/L^2;Z2=6*E*Iy/L^2;Y3=4*E*Iz/L;Z3=4*E*Iy/L;Y4=2*E*Iz/L;
        Z4=2*E*Iy/L;S=G*Ix/L;
    k_loc =[X 0 0 0 0 0 -X 0 0 0 0 0
               0 Y1 0 0 0  Y2 0 -Y1 0 0 0 Y2
               0 0 Z1 0 -Z2 0 0 0 -Z1 0 -Z2 0
               0 0 0 S 0 0 0 0 0 -S 0 0
               0 0 -Z2 0 Z3 0 0 0 Z2 0 Z4 0
               0 Y2 0 0 0 Y3 0 -Y2 0 0 0 Y4
               -X 0 0 0 0 0 X 0 0 0 0 0
               0 -Y1 0 0 0  -Y2 0 Y1 0 0 0 -Y2
               0 0 -Z1 0 Z2 0 0 0 Z1 0 Z2 0
               0 0 0 -S 0 0 0 0 0 S 0 0
               0 0 -Z2 0 Z4 0 0 0 Z2 0 Z3 0
               0 Y2 0 0 0 Y4 0 -Y2 0 0 0 Y3];
       v1=dir/L;
       p1=nodes(elem(i,1),:);
       v2_p= p_aux(i,:)-p1;
       v3=cross(v1,v2_p)/norm(cross(v1,v2_p));
       v2=cross(v3,v1);
       lamda=[v1;v2;v3];
       T=zeros(12);
       T(1:3,1:3)=lamda; T(4:6,4:6)=lamda; T(7:9,7:9)=lamda; T(10:12,10:12)=lamda;
       k_g(ubic,ubic)= k_g(ubic,ubic) + T'*k_loc*T;
    
end

bc=zeros(24,1);
bc(1:6) = ones(6,1); bc(19:24) = ones(6,1);

k_red=k_g(~bc,~bc);
fzas=zeros(24,1);
fzas(8)=-5e3; fzas(10)=-100e3*12; fzas(15)=40e3;
fzas_red=fzas(~bc);

despl=zeros(24,1);
despl(~bc)=k_red^-1*fzas_red;

fzas_calc=k_g*despl;
