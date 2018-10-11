clear;clc;

elementos=load('ElementosT1.txt');
nodos=load('NodosT1.txt');
nodos=nodos*1000;

%% Propiedades del Material
E=30e3;
nu=0.18;
%Plane Strain
C = (E/((1+nu)*(1-2*nu)))*[1-nu    nu      0;
                           nu  1-nu      0;
                            0    0  0.5-nu];
%% Definiciones
nDofNod = 2;                    % grados de libertad por nodo
nel = size(elementos,1);         % elementos
nNod = size(nodos,1);           % nodos
nNodEle = size(elementos,2);     % nodos por elemento
nDofTot = nDofNod*nNod;         % grados de libertad
%% Condiciones de borde
bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
bc(19:24,:) = true;
%% Cargas
R = zeros(nNod,nDofNod);        % Vector de cargas

% Puntos de Gauss
rsInt = 3*ones(1,1);
[wpg, upg, npg] = gauss(rsInt);


%Carga de Compresion superrior
for iele=6:10
    r=0;
    nodesEle=nodos(elementos(iele,:),:);
    for ipg=1:npg
    eta=1;
    ksi=upg(ipg);
    
    N=shapefuns([ksi eta],'Q4');
    
    dN=shapefunsder([ksi eta],'Q4');
    
    jac=dN*nodesEle;
    
    ene(1,1:2:7)=N;
    ene(2,2:2:8)=N; 
    
    q=zeros(4,2);
    q(3:4,2)=-.480;
    q=reshape(q',[8 1]);
    
    I=ene'*ene*q*jac(1,1);
    
    r=r+wpg(ipg)*I;
    
    
    
    end
   R(elementos(iele,:),:)= R(elementos(iele,:),:) + reshape(r,[2 4])';
end
%carga del Agua
carga=@(y) 11/1050000*y -11/20;

for iele=15:5:70
    r=0;
    nodesEle=nodos(elementos(iele,:),:);

   for ipg=1:npg
    ksi=1;
    eta=upg(ipg);
    
    N=shapefuns([ksi eta],'Q4');
    
    dN=shapefunsder([ksi eta],'Q4');
    
    jac=dN*nodesEle;
    
    alturas=nodesEle(2:3,2);
    
    sigmaNodal=[0;carga(alturas);0];
    
    sigma=N*sigmaNodal;
    
    I=[N(2);N(3)]*sigma*jac(2,2);
    
    r=r+I*wpg(ipg);
   
   end
    R(elementos(iele,2:3),1)= R(elementos(iele,2:3),1) + r;
end

%Cargas de volumen

% Puntos de Gauss
rsInt = 3*ones(1,2);
[wpg, upg, npg] = gauss(rsInt);
ro=2e-6; g=-10e3;
for iele=1:nel
   r=0;
   nodesEle=nodos(elementos(iele,:),:);
   for ipg=1:npg
       ksi=upg(ipg,1);
       eta=upg(ipg,2);
    
       N=shapefuns([ksi eta],'Q4');
    
       dN=shapefunsder([ksi eta],'Q4');
       
       jac=dN*nodesEle;
       
       ene(1,1:2:7) = N; ene(2,2:2:8)=N;
       
       Djac = det(jac);
       
       fuerzaNodal=ro*g*[0 1 0 1 0 1 0 1]';
       
       I=ene'*ene*Djac*fuerzaNodal;
       
       
       r=r+I*wpg(ipg);
       
   end
    R(elementos(iele,:),:)= R(elementos(iele,:),:) + reshape(r,[2 4])';
end

%% Matriz de rigidez
K=zeros(nDofTot);
B=zeros(3,8);
for iele=1:nel
    Ke=zeros(nNodEle*nDofNod);
    nodosEle=nodos(elementos(iele,:),:);
    for ipg=1:npg
       ksi=upg(ipg,1);
       eta=upg(ipg,2);
    
       N=shapefuns([ksi eta],'Q4');
    
       dN=shapefunsder([ksi eta],'Q4');
    
       jac=dN*nodesEle;
       
       Djac = det(jac);
       
       dNxy=jac\dN;
       
       B(1,1:2:7)=dNxy(1,:);
       B(2,2:2:8)=dNxy(2,:);
       B(3,1:2:7)=dNxy(2,:); B(3,2:2:8)=dNxy(1,:);
       
       Ke=Ke + B'*C*B*Djac*wpg(ipg);
    end
    eleDofs = node2dof(elementos(iele,:),nDofNod);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
end



% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(R',[],1);

% Solver
Dr = K(isFree,isFree)\Rr(isFree);

% Reconstrucci�n
D = zeros(nDofTot,1);
D(isFree) = D(isFree) + Dr;

D=reshape(D,[2,nNod])';
