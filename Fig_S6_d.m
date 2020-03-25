%% Code for Fig.S6 (d)
%%
clear all
close all

%% Parameters
r=3e3;
Rsp=800e6;
Rsn=900e6;
Rr=45e6;
RM=900e6;
C1=4.5e-12; C2=.65e-6;

W=50e-9; L=8e-6; d= 30e-9;
mun=0.08; lambda=0.03;
cgsm=C1*C2/(C1+C2);
Km=(1/2)*8*(W/L)*mun*cgsm/(L*L);
Vds=0.5;
R0=1e3;

VD=0.3;
VT0=0.35;
f=0.1;
V1=9; W1=0.25; T1=0.5; N1=1800;


%% Time
dt=0.01;
t=[0:dt:N1*T1-dt];

%% Signal inputs
tadv=(T1-W1)/2;
lhl=1+sin(2*pi/T1*(-T1/4+tadv));
svis=(1+sin(2*pi/T1*(t(1:T1/dt)-T1/4+tadv)));
svip=V1*rectangularPulse(lhl, 2.001, svis);
vip=repmat(svip, 1,N1);
vi=vip;

%% Current calculation

vM(1)=0;
vg(1)=0;
vA(1)=0;
ID(1)=0;
Krhi=abs(r/(R0*abs(V1(1))));
for n=2:length(t)
    if vi(n)~=0
        Krh=abs(r/(R0*abs(vi(n))));
        Krhi=Krh;
    else
        Krh=Krhi;
    end
    if vi(n)>0
        vA(n)=vi(n)/(1+Krh*(vi(n)-vg(n-1)));
    else
        if vi(n)<0
            vA(n)=vi(n)/(1+Krh*(-vi(n)+vg(n-1)));
        else
            vA(n)=vi(n)/(1+Krh*abs(vg(n-1)));
        end
    end
    
    if (vA(n)>VD+vg(n-1)) && (vi(n)>=0)
        vg(n)=vg(n-1)*exp(-dt*(1/Rr+1/Rsp)/C1)+vM(n-1)*(exp(-dt/(RM*C2))-exp(-dt*(1/Rr+1/Rsp)/C1))+...
            Rr/(Rsp+Rr)*(vA(n)-VD-vM(n-1))*(1-exp(-dt*(1/Rr+1/Rsp)/C1))+...
            RM/Rsp*(vA(n)-VD-vg(n-1))*(1-exp(-dt/(RM*C2)));
        
        vM(n)=vM(n-1)*exp(-dt/(RM*C2))+ RM/Rsp*(vA(n)-VD-vg(n))*(1-exp(-dt/(RM*C2)));
    else
        if (vA(n)+VD<vg(n-1)) && (vi(n)<=0)
            
            vg(n)=vg(n-1)*exp(-dt*(1/Rr+1/Rsn)/C1)+vM(n-1)*(exp(-dt/(RM*C2))-exp(-dt*(1/Rr+1/Rsn)/C1))+...
                Rr/(Rsn+Rr)*(vA(n)+VD-vM(n-1))*(1-exp(-dt*(1/Rr+1/Rsn)/C1))+...
                RM/Rsn*(vA(n)+VD-vg(n-1))*(1-exp(-dt/(RM*C2)));
            
            vM(n)=vM(n-1)*exp(-dt/(RM*C2))+ RM/Rsn*(vA(n)+VD-vg(n))*(1-exp(-dt/(RM*C2)));
        else
            vg(n)=vM(n-1)*exp(-dt/(RM*C2))+(vg(n-1)-vM(n-1))*exp(-dt/(Rr*C1));
            vM(n)=vM(n-1)*exp(-dt/(RM*C2));
        end
    end
    vTH(n)=VT0*exp(-f*abs(vM(n))/VT0);
    
    if (vg(n)-vTH(n))>Vds
        ID(n)=Km*(2*(vg(n)-vTH(n))-Vds)*Vds;
    else
        if vg(n)>vTH(n)
            ID(n)=Km*(vg(n)-vTH(n))*(vg(n)-vTH(n));
        else
            ID(n)=0;
        end
    end
end

figure; plot(t,ID);
