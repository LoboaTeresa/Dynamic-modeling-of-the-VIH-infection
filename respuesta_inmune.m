clear all 
close all
clc

Vo=10^-6;
To=10; %porcentaje de celulas CD4+ T en la sangre periferica;
Io=0;

%Parametro paciente 2:

d=0.020
lambda=d*To;
c=3;
k=0.36*10^-3;%0.65*10^-3;
delta0=0.80;
p=1800;

beta=0.090;
t2=65;
t1=60;

%Parametros fijos:
kp=1+beta*10^5;
inT1=2.5;
inT2=15;


to=20;
tn=400;
N=10000;
%Lo incluimos en la representacion
h = @(T,I,V,t)( lambda - k.*V.*T - d.*T );
g = @(T,I,V,t)( k.*V.*T - ( delta0 + (    (beta ./ ( 1 + kp.*exp( -(t-t1)./inT1 )) )   -   (beta ./ ( 1 + kp.*exp(-(t-t2)./inT2) )  ) ).*V).*I );
f = @(T,I,V,t)( p.*I - c.*V );

[Tr,Ir,Vr,tr]=lorungek4(h,g,f,to,To,Io,Vo,tn,N);

h = @(T,I,V,t)( lambda - k.*V.*T - d.*T );
g = @(T,I,V,t)( k.*V.*T - delta0.*I );
f = @(T,I,V,t)( p.*I - c.*V );

[Ti,Ii,Vi,ti]=lorungek4(h,g,f,to,To,Io,Vo,tn,N);


figure(1)
plot(tr,log(Vr)+14); title('Modelo de infecci?n v?rica con respuesta inmune'); xlabel('Tiempo / (Dias)'); ylabel('VIH RNA ml^-^1')
hold on
plot(ti,log(Vi)+14,'r'); legend('Modelo con respuesta inmune','Modelo sin respuesta inmune')
plot(35,8,'rx','MarkerSize',12)
plot(40,14,'xr','MarkerSize',12)
plot(45,18,'xr','MarkerSize',12)
plot(90,16,'xr','MarkerSize',12)
plot(120,16,'xr','MarkerSize',12)
plot(170,16,'xr','MarkerSize',12)
plot(250,18,'xr','MarkerSize',12)
plot(350,16,'xr','MarkerSize',12)
hold off

