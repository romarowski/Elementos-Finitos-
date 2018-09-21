E=210e3; nu=.3; G= E/(2*(1+nu)); h=150; b=20; d=50;
Iz=b*h^3/12; Iy=b^3*h/12; K=.229*h*b^3; A=h*b;
I=pi*d^4/64; Kp=pi*d^4/32; Ap=pi*d^2/4;

dofs=42;

nodeDofs=1:42; nodeDofs=reshape(nodeDofs,6,7); nodeDofs=nodeDofs';

elem=[1 2;3 4;5 6;2 4;2 6;2 7;4 6;4 7;6 7];
nodes=[0 0 0;0  0 1500;500 866.03 0; 500 866.03 1500;-500 866.03 0;-500 866.03 1500;0 577.35 2316.5];
p_3=[-10 0 0;0 577.35 0;0 577.35 0];
p_p=[0 577.35 1500];
k_g=zeros(dofs);
m=1;
for i=1:9
    if i>3
    Iz=I; Iy=I; K=Kp; A=Ap;
    p_3=p_p;
    end
    ubic=[nodeDofs(elem(i,1),:) nodeDofs(elem(i,2),:)];
    L=norm(nodes(elem(i,2),:)-nodes(elem(i,1),:));
    dir=nodes(elem(i,2),:)-nodes(elem(i,1),:);
    
X=A*E/L;Y1=12*E*Iz/L^3;Z1=12*E*Iy/L^3;Y2=6*E*Iz/L^2;Z2=6*E*Iy/L^2;Y3=4*E*Iz/L;Z3=4*E*Iy/L;Y4=2*E*Iz/L;
Z4=2*E*Iy/L;S=G*K/L;

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
    if i>3
        v2_p=p_3 -p1;
    else
        v2_p= p_3(m,:) -p1;
        m=m+1;
    end
    v3=cross(v1,v2_p)/norm(cross(v1,v2_p));
    v2=cross(v3,v1);
    lamda=[v1;v2;v3];
    T=zeros(12);
    T(1:3,1:3)=lamda; T(4:6,4:6)=lamda; T(7:9,7:9)=lamda; T(10:12,10:12)=lamda;
    k_g(ubic,ubic)= T'*k_loc*T +k_g(ubic,ubic);
end

bc=zeros(42,1);
bc(1:6)=ones(6,1); bc(25:30)=ones(6,1); bc(13:18)=ones(6,1);

fzas=zeros(42,1); fzas(37)=1.5e3; fzas(34)=-500; fzas(35)=866.03; fzas(10)=-500; fzas(11)=866.03;
fzas_red=fzas(~bc);
k_red=k_g(~bc,~bc);
despl=zeros(42,1);

despl(~bc)=k_red^(-1)*fzas_red;

i=5;

ubic=[nodeDofs(elem(i,1),:) nodeDofs(elem(i,2),:)];
L=norm(nodes(elem(i,2),:)-nodes(elem(i,1),:));
dir=nodes(elem(i,2),:)-nodes(elem(i,1),:);
v1=dir/L;
p1=nodes(elem(i,1),:);
v2_p=p_3 -p1;
v3=cross(v1,v2_p)/norm(cross(v1,v2_p));
v2=cross(v3,v1);
lamda=[v1;v2;v3];
T=zeros(12);
T(1:3,1:3)=lamda; T(4:6,4:6)=lamda; T(7:9,7:9)=lamda; T(10:12,10:12)=lamda;
despl_loc=T*despl(ubic);
B_bend=[-6/L^2+12*L/2/L/L^3 -4/L+6*L/2/L^2 6/L^2-12*L/2/L/L^3 -2/L+6*L/2/L^2];
B_axial=[-1/L 1/L];
sigma_b_z= E*d/2*B_bend*despl_loc([2 6 8 12]);
sigma_b_y= E*d/2*B_bend*despl_loc([3 5 9 11]);
sigma_axial=E*B_axial*despl_loc([1 7]);
sigma_tot=sigma_axial+sigma_b_y+sigma_b_z;
B_shear=[12/L^3 6/L^2 -12/L^3 6/L^2];
tau_shear_xy=2/3*d^3/8*E/d*B_shear*despl_loc([2 6 9 12]);
tau_shear_xz=2/3*d^3/8*E/d*B_shear*despl_loc([3 5 9 11]);
tau=G*d/2*B_axial*despl_loc([4 10]);

