clear all 
close all
clc

%%Plasma Virus en Infeccion primaria de VIH

%Condiciones iniciales:

Vo=10^-6;
To=1; %porcentaje;
Io=0;

%Par?metros(media):

d=0.01;
lambda=d*To;
c=3;
k=0.65;
delta=0.39;
p=850;

%Ecuaciones diferenciales:

h = @(T,I,V,t)( lambda - k.*V.*T - d.*T );
g = @(T,I,V,t)( k.*V.*T - delta.*I );
f = @(T,I,V,t)( p.*I - c.*V );

[T,I,V,t]=rungek23(h,g,f,to,To,Io,Vo,tn,N);
