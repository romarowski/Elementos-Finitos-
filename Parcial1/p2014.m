E=70e3; b=50; h=100; Iz=h^3*b/12; A=b*h;

nodeDofs=[1 2 3; 4:6;7:9];
nodes=[0 0;-0.5 0;-0.75 0]; nodes=nodes*10^3;
elem=[1 2;2 3];

k_g=zeros(9);
for i=1:2
    ubic=[nodeDofs(elem(i,1),:) nodeDofs(elem(i,2),:)];
    Le=norm( nodes(elem(i,2),:) - nodes(elem(i,1),:));
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
    dts= (nodes(elem(i,2),:) - nodes(elem(i,1),:))/Le;
    T=zeros(6); T(1,1:2)=dts; T(4,4:5)=dts; T(2,1:2)=[-1*dts(2) dts(1)]; T(5,4:5)=[-1*dts(2) dts(1)];
    T(3,3)=1; T(6,6)=1;
    k_g(ubic,ubic)= T'*Kel*T + k_g(ubic,ubic);
    
end

bc=zeros(9,1); bc(1:2)=[1 1]; bc(4:5)=[1 1]; 

k_red=k_g(~bc,~bc);

fzas=zeros(9,1);
fzas(8)=-.125e3; fzas(9)=-5.2e3; fzas(5)=(-.125-.25)*1e3; fzas(6)=(5.2e-3-.0208)*1e6; fzas(2)=-.25e3; 
fzas(3)=0.0208e6 -1;

fzas_red=fzas(~bc);
despl=zeros(9,1);

despl(~bc)= k_red^-1*fzas_red;
despl_otra_mitad= despl; 
despl_otra_mitad(8)=despl_otra_mitad(8)*-1;

x=0.5e3;
L=0.5e3;
N=[1-3*x^2/L^2+2*x^3/L^3;x-2*x^2/L+x^3/L^2; 3*x^2/L^2-2*x^3/L^3;-x^2/L+x^3/L^2];
v1= N'*despl([2 3 5 6]);
L=.25e3;
N=[1-3*x^2/L^2+2*x^3/L^3;x-2*x^2/L+x^3/L^2; 3*x^2/L^2-2*x^3/L^3;-x^2/L+x^3/L^2];
v2= N'*despl([5 6 8 9]);
if x>0.5
    v=v1;
else
    v=v2;
end

L=.25e3; x=125;
B = [-6/L^2 + 12*x/L^3 ; -4/L+6*x/L^2; 6/L^2 - 12*x/L^3; -2/L+6*x/L^2];
B_ax= [-1/L 1/L];
sig_ax=E*B_ax*[despl(4) despl(7)]';
sig_bend= E*50*B'*[despl(8) despl(9) despl(5) despl(6)]';

