I = 83.33e5; Aviga=10e3; Abarra = 1.96e3; E=210e3;

elem=[1 3
      5 6
      3 4
      6 7
      3 6
      7 4
      2 3];
k_g = zeros(21);
  %sumo las vigas
for i=1:2
    ubic_global = [3*elem(i,1)-2:3*elem(i,1)  3*elem(i,2)-2:3*elem(i,2) ];
    k_g(ubic_global,ubic_global) = elemviga(E,I,Aviga,400,1,deg2rad(36.87)) + k_g(ubic_global,ubic_global); 
 end
for i=3:4
    ubic_global = [3*elem(i,1)-2:3*elem(i,1)  3*elem(i,2)-2:3*elem(i,2) ];
    k_g(ubic_global,ubic_global) = elemviga(E,I,Aviga,500,1,deg2rad(36.87)) + k_g(ubic_global,ubic_global); 
end
%sumo las barras
ubic_global = [3*elem(5,1)-2:3*elem(5,1)-1  3*elem(5,2)-2:3*elem(5,2)-1 ];
k_g(ubic_global,ubic_global) = elembarra(1,E,Abarra,240,deg2rad(90)) + k_g(ubic_global,ubic_global);
ubic_global = [3*elem(6,1)-2:3*elem(6,1)-1  3*elem(6,2)-2:3*elem(6,2)-1 ];
k_g(ubic_global,ubic_global) = elembarra(1,15*E,Abarra,240,deg2rad(90)) + k_g(ubic_global,ubic_global);
ubic_global = [3*elem(7,1)-2:3*elem(7,1)-1  3*elem(7,2)-2:3*elem(7,2)-1 ];
k_g(ubic_global,ubic_global) = elembarra(1,E,Abarra,424.3,deg2rad(45)) + k_g(ubic_global,ubic_global);
 
bc = zeros(21,1);
bc(1:2) = [1 1]; bc(4:5) = [1 1]; bc(13:14) = [1 1];
k_red = k_g(~bc,~bc);
 
fzas = zeros(21,1);
fzas(11)=-4e3; fzas(12)= -112.5e4; fzas(20) = -4e3; fzas(21) = -112.5e4;
fzas_red = fzas(~bc);
  
despl= k_red*fzas_red;
  
  
  