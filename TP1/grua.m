E=210e3; diam_barra=5; lado_viga=10; 
area_barra=pi*diam_barra^2/4; area_viga=lado_viga^2; I_beam=lado_viga^4/12;
dofs=23; 

youngs=[E E E 100*E];
nodeDofs =[1 2 0;3 4 5;6 7 8;9 10 11;9 10 12;13:15;16:18;19 20 0;21:23];
elem=[1 3;2 3;3 4;5 6;6 7;7 8;6 9;5 9;8 9];
nodes=[0 0;0 24;-30 30;-51 42;-51 42;-46 57;-46 65;-66 65;-66 49];
isBarra=[1 0 0 0 0 1 1 0 1]; isBarra= logical(isBarra);

k_g=zeros(dofs); j=1;
for i=1:length(elem)
    if isBarra(i)
        ubic= [nodeDofs(elem(i,1),1:2) nodeDofs(elem(i,2),1:2)];
        L = norm( nodes(elem(i,2),:) - nodes(elem(i,1),:));
        k_loc = area_barra*youngs(j)/L * [1 -1;...
                             -1 1];          
        dts = (nodes(elem(i,2),:)-nodes(elem(i,1),:))/L ;
        T=[dts 0 0; 0 0 dts];
        k_g(ubic,ubic) = T'*k_loc*T + k_g(ubic,ubic);
        j=j+1;
    else
        ubic = [nodeDofs(elem(i,1),:),nodeDofs(elem(i,2),:)];
        L = norm( nodes(elem(i,2),:) - nodes(elem(i,1),:));
        X = area_viga*E/L; Y_1=12*E*I_beam/L^3; Y_2=6*E*I_beam/L^2; Y_3=4*E*I_beam/L; Y_4 = 2*E*I_beam/L;  
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
bc(1:2)=ones(2,1); bc(3:4)=ones(2,1);

k_red= k_g(~bc,~bc);

fzas=zeros(dofs,1);
fzas(19)=3000; fzas(20)=-1600; fzas(21)=-3000; fzas(22)=-1600;
fzas_red=fzas(~bc);

despl=zeros(dofs,1);

despl(~bc) =k_red^(-1)*fzas_red; 



