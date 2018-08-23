E=200e3; L=1e3;  diam=40; a=40;
I_bar = pi*diam^4/64; I_beam = a^4/12; A_bar = pi*diam^4/4; A_beam=a^2;

elem =[1 2; 3 8; 4 9; 4 5; 6 10; 3 7];

nodes= [0 0; 0.71*L .71*L ; 1.71*L .71*L ;1.353*L .353*L ; L 0; 2*L 0; 2.71*L .71 *L];

BARRA = true;
esbarra=[BARRA 0 0 0 BARRA 0];
nodeDofs=[1 2 0; 3 4 5;6 7 8;12 13 14 ;15 16 17; 18 19 0; 9 10 11;3 4 21 ;6 7 20;12 13 21];







