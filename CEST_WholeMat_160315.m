% CEST experiment for exchange between two states, J coupling, Full simulation.

% Necessary Definitions
GammaB1=10;
omegaIp=GammaB1*pi*2*exp(1i*(pi/2));
omegaIm=GammaB1*pi*2*exp(-1i*(pi/2));
J=0;

omegaSp=0;
omegaSm=0;

NBF=86.16;  % N base frequency in MHz

StateBshift=5;

R1Sa=-2;
R1Sb=-2;
R1Ia=-1;
R1Ib=-1;

R2Sa=-12;
R2Sb=-12;
R2Ia=-18;
R2Ib=-18;
Rmq=-30;

sigma=0;
delS=0;
delI=0;
etaS=0;
mu_mq=0;
CESTperiod=0.5;

pb=0.1;
pa=1-pb;
kex=250;

kab=pb*kex;
kba=pa*kex;


% Initial State Vector
I0a=[0           0               pa/2                       0        0            pa/2      0       0       0       0       1/2       0        0          0       0]';
I0b=[0           0               pb/2                       0        0            pb/2      0       0       0       0       1/2       0        0          0       0]';
I0=[I0a;
    I0b];
% offset of CEST RF
Ia_ppmOffsetRanges=-15:0.1:15;
clear CEST_DAT;
CEST_DAT=zeros(length(Ia_ppmOffsetRanges), 6);

% offset of CEST RF
Ib_ppmOffset=Ia_ppmOffsetRanges-StateBshift;

OmegaSa=50;
OmegaSb=50;
zeromat15=zeros(15, 15);
Onemat15=eye(15);

for counter1=1:length(Ia_ppmOffsetRanges)

OmegaIa=Ia_ppmOffsetRanges(counter1)*NBF*2*pi;
OmegaIb=Ib_ppmOffset(counter1)*NBF*2*pi;

d1=R2Ia+0.5*R1Sa-1i*OmegaIa-1i*pi*J;
d2=R2Ia+0.5*R1Sa+1i*OmegaIa+1i*pi*J;
d3=R1Ia+0.5*R1Sa;
d4=R2Ia+0.5*R1Sa-1i*OmegaIa+1i*pi*J;
d5=R2Ia+0.5*R1Sa+1i*OmegaIa-1i*pi*J;
d6=R1Ia+0.5*R1Sa;
d7=R2Sa+0.5*R1Ia+etaS-1i*(OmegaSa+pi*J);
d8=R2Sa+0.5*R1Ia+etaS+1i*(OmegaSa+pi*J);
d9=R2Sa+0.5*R1Ia+etaS-1i*(OmegaSa-pi*J);
d10=R2Sa+0.5*R1Ia+etaS+1i*(OmegaSa-pi*J);
d11=R1Sa;
d12=Rmq-mu_mq-1i*(OmegaIa-OmegaSa);
d13=Rmq-mu_mq+1i*(OmegaIa-OmegaSa);
d14=Rmq+mu_mq-1i*(OmegaIa+OmegaSa);
d15=Rmq+mu_mq+1i*(OmegaIa+OmegaSa);

% General relaxation matrix
Ra=[d1           0                   1i*omegaIp        -R1Sa/2      0           0           0           0           0           0           0           -1i*omegaSp/2     0            1i*omegaSm/2      0;
    0           d2                  -1i*omegaIm       0           -R1Sa/2      0           0           0           0           0           0           0     1i*omegaSm/2   0               -1i*omegaSp/2;
    1i*omegaIm/2  -1i*omegaIp/2         d3              0           0           -R1Sa/2      1i*omegaSm/4  -1i*omegaSp/4 -1i*omegaSm/4 1i*omegaSp/4  (sigma+delS)/2  0           0            0               0;
    -R1Sa/2      0                   0               d4          0           1i*omegaIp    0           0           0           0           0           1i*omegaSp/2      0            -1i*omegaSm/2     0;
    0           -R1Sa/2              0               0           d5          -1i*omegaIm   0           0           0           0           0           0        -1i*omegaSm/2  0               1i*omegaSp/2;
    0           0                   -R1Sa/2          1i*omegaIm/2  -1i*omegaIp/2 d6   -1i*omegaSm/4 1i*omegaSp/4  1i*omegaSm/4 1i*omegaSp/4 (sigma-delS)/2  0       0             0                0;
    0           0                   1i*omegaSp/2     0           0           -1i*omegaSp/2 d7         0           -R1Ia/2       0          1i*omegaSp/2  0       -1i*omegaIp/2 -1i*omegaIm/2      0;
    0           0                   -1i*omegaSm/2    0           0           1i*omegaSm/2  0          d8          0           -R1Ia/2      -1i*omegaSm/2 1i*omegaIm/2 0          0               -1i*omegaIp/2;
    0           0                   -1i*omegaSp/2    0           0           1i*omegaSp/2  -R1Ia/2     0           d9          0           1i*omegaSp/2  0          1i*omegaIp/2 -1i*omegaIm/2    0;
    0           0                   1i*omegaSm/2      0           0          -1i*omegaSm/2 0          -R1Ia/2      0           d10         -1i*omegaSm/2 -1i*omegaIm/2 0         0               1i*omegaIp/2;
    0           0                   sigma+delS      0           0           sigma-delS  1i*omegaSm/2 -1i*omegaSp/2  1i*omegaSm/2 -1i*omegaSp/2  d11      0          0           0               0;
    -1i*omegaSm/2  0                 0           1i*omegaSm/2     0           0           0           1i*omegaIp/2 0           -1i*omegaIp/2 0          d12         0           0               0;
    0           1i*omegaSp/2         0           0              -1i*omegaSp/2 0           -1i*omegaIm/2    0       1i*omegaIm/2 0           0           0           d13         0               0;
    1i*omegaSp/2  0                 0           -1i*omegaSp/2     0           0           1i*omegaIp/2    0       -1i*omegaIp/2 0           -1i*omegaIp/2 0          0           d14             0;
    0           -1i*omegaSm/2        0           0           1i*omegaIm/2     0           0              -1i*omegaIm/2 0       1i*omegaIm/2     0        0          0           0               d15];


d1=R2Ib+0.5*R1Sb-1i*OmegaIb-1i*pi*J;
d2=R2Ib+0.5*R1Sb+1i*OmegaIb+1i*pi*J;
d3=R1Ib+0.5*R1Sb;
d4=R2Ib+0.5*R1Sb-1i*OmegaIb+1i*pi*J;
d5=R2Ib+0.5*R1Sb+1i*OmegaIb-1i*pi*J;
d6=R1Ib+0.5*R1Sb;
d7=R2Sb+0.5*R1Ib+etaS-1i*(OmegaSb+pi*J);
d8=R2Sb+0.5*R1Ib+etaS+1i*(OmegaSb+pi*J);
d9=R2Sb+0.5*R1Ib+etaS-1i*(OmegaSb-pi*J);
d10=R2Sb+0.5*R1Ib+etaS+1i*(OmegaSb-pi*J);
d11=R1Sb;
d12=Rmq-mu_mq-1i*(OmegaIb-OmegaSb);
d13=Rmq-mu_mq+1i*(OmegaIb-OmegaSb);
d14=Rmq+mu_mq-1i*(OmegaIb+OmegaSb);
d15=Rmq+mu_mq+1i*(OmegaIb+OmegaSb);
        
Rb=[d1           0                   1i*omegaIp        -R1Sb/2      0           0           0           0           0           0           0           -1i*omegaSp/2     0            1i*omegaSm/2      0;
    0           d2                  -1i*omegaIm       0           -R1Sb/2      0           0           0           0           0           0           0     1i*omegaSm/2   0               -1i*omegaSp/2;
    1i*omegaIm/2  -1i*omegaIp/2         d3              0           0           -R1Sb/2      1i*omegaSm/4  -1i*omegaSp/4 -1i*omegaSm/4 1i*omegaSp/4  (sigma+delS)/2  0           0            0               0;
    -R1Sb/2      0                   0               d4          0           1i*omegaIp    0           0           0           0           0           1i*omegaSp/2      0            -1i*omegaSm/2     0;
    0           -R1Sb/2              0               0           d5          -1i*omegaIm   0           0           0           0           0           0        -1i*omegaSm/2  0               1i*omegaSp/2;
    0           0                   -R1Sb/2          1i*omegaIm/2  -1i*omegaIp/2 d6   -1i*omegaSm/4 1i*omegaSp/4  1i*omegaSm/4 1i*omegaSp/4 (sigma-delS)/2  0       0             0                0;
    0           0                   1i*omegaSp/2     0           0           -1i*omegaSp/2 d7         0           -R1Ib/2       0          1i*omegaSp/2  0       -1i*omegaIp/2 -1i*omegaIm/2      0;
    0           0                   -1i*omegaSm/2    0           0           1i*omegaSm/2  0          d8          0           -R1Ib/2      -1i*omegaSm/2 1i*omegaIm/2 0          0               -1i*omegaIp/2;
    0           0                   -1i*omegaSp/2    0           0           1i*omegaSp/2  -R1Ib/2     0           d9          0           1i*omegaSp/2  0          1i*omegaIp/2 -1i*omegaIm/2    0;
    0           0                   1i*omegaSm/2      0           0          -1i*omegaSm/2 0          -R1Ib/2      0           d10         -1i*omegaSm/2 -1i*omegaIm/2 0         0               1i*omegaIp/2;
    0           0                   sigma+delS      0           0           sigma-delS  1i*omegaSm/2 -1i*omegaSp/2  1i*omegaSm/2 -1i*omegaSp/2  d11      0          0           0               0;
    -1i*omegaSm/2  0                 0           1i*omegaSm/2     0           0           0           1i*omegaIp/2 0           -1i*omegaIp/2 0          d12         0           0               0;
    0           1i*omegaSp/2         0           0              -1i*omegaSp/2 0           -1i*omegaIm/2    0       1i*omegaIm/2 0           0           0           d13         0               0;
    1i*omegaSp/2  0                 0           -1i*omegaSp/2     0           0           1i*omegaIp/2    0       -1i*omegaIp/2 0           -1i*omegaIp/2 0          0           d14             0;
    0           -1i*omegaSm/2        0           0           1i*omegaIm/2     0           0              -1i*omegaIm/2 0       1i*omegaIm/2     0        0          0           0               d15];
     
            
        
        
        
        
Rab=[Ra  zeromat15;
    zeromat15 Rb];
exchangemat=[-kab   kba;
              kab   -kba];
          exchangemat=kron(exchangemat,Onemat15);
          
          Rab=Rab+exchangemat;
        
        
        
        
        M=expm(Rab*CESTperiod)*I0;
        
        CEST_DAT(counter1,1)=OmegaIa/(NBF*2*pi);
        CEST_DAT(counter1,2)=M(3);
        CEST_DAT(counter1,3)=M(6);
        CEST_DAT(counter1,4)=M(18);
        CEST_DAT(counter1,5)=M(21);
        CEST_DAT(counter1,6)=M(3)+M(6);

end;


figure();
hold all;
plot(CEST_DAT(:,1), CEST_DAT(:,6));