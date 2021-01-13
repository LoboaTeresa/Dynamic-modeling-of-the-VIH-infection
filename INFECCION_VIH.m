clear all 
close all
clc

%%Plasma Virus en Infeccion primaria de VIH

%Condiciones iniciales:

Vo=10^-6;
To=10; %porcentaje de celulas CD4+ T en la sangre periferica;
Io=0;

%Par?metros(media):

d=0.012;%0.01;
lambda=d*To;
c=3;
k=0.75*10^-3;%0.65*10^-3;
delta=0.39;
p=790;

%Periodo de representaci?n:

to=20;
tn=400;
N=10000;

%Ecuaciones diferenciales:

h = @(T,I,V,t)( lambda - k.*V.*T - d.*T );
g = @(T,I,V,t)( k.*V.*T - delta.*I );
f = @(T,I,V,t)( p.*I - c.*V );

[T,I,V,t]=lorungek4(h,g,f,to,To,Io,Vo,tn,N);

figure(1)

plot(t,log(V)+14); title('Infecci?n primaria del virus del VIH');xlabel('Tiempo / (d?as)') ;ylabel('HIV-1 RNA / (ml^-^1)');

hold on
A=ones(1,36).*( (p*lambda)/(delta*c) -(d/k) );
A=A-46.9;
plot(t(to:280:end),A,'.-.');legend('Virus VIH','Valor estable virus');
plot(35,8,'rx','MarkerSize',12)
plot(40,14,'xr','MarkerSize',12)
plot(45,18,'xr','MarkerSize',12)
plot(90,16,'xr','MarkerSize',12)
plot(120,16,'xr','MarkerSize',12)
plot(170,16,'xr','MarkerSize',12)
plot(250,18,'xr','MarkerSize',12)
plot(350,16,'xr','MarkerSize',12)
hold off
%%
R= (k*p*lambda)/(c*d*delta); %ratio de celulas infectadas que genera una celula infectada antes de morir
Tdoble=log(2)/( delta*(R-1) ); %tiempo en doblarse la cantidad de virus

%Para periodos mayores datos por debajo de niveles estables:
%Cambia delta;

%Par?metros paciente 2:

beta=0.090;
t2=65;
t1=60;
delta0=0.80;

%Parametros fijos:
k=1+beta*10^5;
inT1=2.5;
inT2=15;

%funcion que varia beta a partir de t1


T1=linspace(to,t1,85);
T2=linspace(t1,tn,N-86);

fun=@(T2)( ( beta ./ ( 1 + k.*exp(-(T2-t1)./inT1) ) ) - ( beta ./ ( 1 + k.*exp(-(T2-t2)./inT2) ) ) );

figure(2)
plot(T1,delta0,'r'); title('Evoluci?n de \delta por la respuesta inmune');xlabel('Tiempo / (d?as)'); ylabel('\delta')
hold on

plot(T2,delta0+fun(T2).*log(V(87:end)))
hold off




%% Representamos segundo modelo de respuesta inmune
d=0.012;%0.01;
lambda=d*To;
c=3;
k=0.75*10^-3;%0.65*10^-3;
delta=0.39;
p=790;

%Periodo de representaci?n:

to=20;
tn=500;
N=10000;

%Ecuaciones diferenciales:

h = @(T,I,V,t)( lambda - k.*V.*T - d.*T );
g = @(T,I,V,t)( k.*V.*T - delta.*I );
f = @(T,I,V,t)( p.*I - c.*V );

[T,I,V,t]=lorungek4(h,g,f,to,To,Io,Vo,tn,N);
%Cambiamos p:
p=790*0.6
h = @(T,I,V,t)( lambda - k.*V.*T - d.*T );
g = @(T,I,V,t)( k.*V.*T - delta.*I );
f = @(T,I,V,t)( p.*I - c.*V );

[Ti,Ii,Vi,ti]=lorungek4(h,g,f,to,To,Io,Vo,tn,N);

figure(4)

plot(t,log(V)+14); title('Respuesta inmune inhibiendo la producci?n de virus');xlabel('Tiempo / (dias)') ;ylabel('HIV-1 RNA / (ml^-^1)');

hold on
plot(ti,log(Vi)+14,'-r'); legend('Modelo sin respuesta inmune', 'Modelo con respuesta inmune')

plot(40,10,'rx','MarkerSize',12)
plot(40,15,'xr','MarkerSize',12)
plot(50,20,'xr','MarkerSize',12)
plot(90,15,'xr','MarkerSize',12)
plot(150,16,'xr','MarkerSize',12)
plot(170,17,'xr','MarkerSize',12)
plot(250,15,'xr','MarkerSize',12)
plot(350,15,'xr','MarkerSize',12)
plot(450,13,'xr','MarkerSize',12)
hold off