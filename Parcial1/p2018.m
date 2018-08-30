E=210e3; P=1.5e3; M_x=-499.99; M_y = 866.0266;
G= E/2/(1+.3);
diam=50; a=20; b=150; 

nodeDofs=1:42; nodeDofs=reshape(nodeDofs,6,7); nodeDofs=nodeDofs';

elem=[1 7;3 4;2 6;4 7;5 7;4 5;4 6;5 6;6 7];
nodes=[0 0 0;-500 866.03 0;500 866.03 0;500 866.03 1500;0 577.35 2316.5;-500 866.03 1500;0 0 1500];

isCirc=logical([0 0 0 ones(1,6)]');
dofs=42;
vecdep3=[0 10 0; 0 866.03 0;0 866.03 0];
k_g=zeros(dofs);
j=1;

for i=1:length(elem)
    if isCirc(i)
        A=pi*diam^2/4; Iy=pi*diam^4/64; Iz=Iy; Ip=Iy*2;
        ubic=[nodeDofs(elem(i,1),:) nodeDofs(elem(i,2),:)];
        dir=nodes(elem(i,2),:) - nodes(elem(i,1),:);
        L=norm(nodes(elem(i,2),:) - nodes(elem(i,1),:));
        X=A*E/L;Y1=12*E*Iz/L^3;Z1=12*E*Iy/L^3;Y2=6*E*Iz/L^2;Z2=6*E*Iy/L^2;Y3=4*E*Iz/L;Z3=4*E*Iy/L;Y4=2*E*Iz/L;
        Z4=2*E*Iy/L;S=G*Ip/L;
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
           p3=[0 0 0];  v2_prima = p3 - p1;
           v3=cross(v1,v2_prima)/norm(cross(v1,v2_prima));
           v2=cross(v3,v1);
           lamda =[v1;v2;v3];
           T=zeros(12); T(1:3,1:3)=lamda; T(4:6,4:6)=lamda; T(7:9,7:9)=lamda; T(10:12,10:12)=lamda;
           k_g(ubic,ubic) = k_g(ubic,ubic)+ T'*k_loc*T;
           %isnan(T'*k_loc*T)
           %i

    else
        A=a*b; Iz=a*b^3/12; Iy=a^3*b/12; Ip=a*b/12*(a^2+b^2);
        ubic=[nodeDofs(elem(i,1),:) nodeDofs(elem(i,2),:)];
        dir=nodes(elem(i,2),:) - nodes(elem(i,1),:);
        L=norm(nodes(elem(i,2),:) - nodes(elem(i,1),:));
        X=A*E/L;Y1=12*E*Iz/L^3;Z1=12*E*Iy/L^3;Y2=6*E*Iz/L^2;Z2=6*E*Iy/L^2;Y3=4*E*Iz/L;Z3=4*E*Iy/L;Y4=2*E*Iz/L;
        Z4=2*E*Iy/L;S=G*Ip/L;
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
           p3=vecdep3(j,:);
           v2_prima = p3 - p1;
           v3=cross(v1,v2_prima)/norm(cross(v1,v2_prima));
           v2=cross(v3,v1);
           lamda =[v1;v2;v3];
           T=zeros(12); T(1:3,1:3)=lamda; T(4:6,4:6)=lamda; T(7:9,7:9)=lamda; T(10:12,10:12)=lamda;
           k_g(ubic,ubic) = k_g(ubic,ubic)+ T'*k_loc*T;
           j=j+1;        
    end
    
end

bc=zeros(dofs,1);
bc(1:18)= ones(18,1);

k_red=k_g(~bc,~bc);
fzas=zeros(dofs,1);
fzas(25)=1.5e3; fzas(34:35)=[-499.99 866.0266]; fzas(40:41)=[-499.99 866.0266];
fzas_red=fzas(~bc);

despl=zeros(dofs,1);

despl(~bc)=k_red^(-1)*fzas_red;
