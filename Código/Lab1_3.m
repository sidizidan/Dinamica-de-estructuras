clc;close all;clear 

%Cargar archivos
g11 = load('g21.txt'); %Registro Desplazamiento incial
g12 = load('g22.txt'); %Registro Velocidad inicial
g13 = load('g23.txt'); %Registro Amortiguamiento adicional
g14 = load('g24.txt'); %Registro Resonancia 

%Constantes sensores
c_respuesta = 1.028; 
c_forzante = 0.986;

%Vectores de tiempo de cada registro
t1 = g11(:,1);
t2 = g12(:,1);
t3 = g13(:,1);
t4 = g14(:,1);

%Respuesta aceleración sistema [g]
Respuesta_A_1 = g11(:,2)./c_respuesta;
Respuesta_A_2 = g12(:,2)./c_respuesta;
Respuesta_A_3 = g13(:,2)./c_respuesta;
Respuesta_A_4 = g14(:,2)./c_respuesta;

%Gráfico respuestas aceleración
figure()
nexttile
plot(t1,Respuesta_A_1)
xlabel('tiempo [s]')
ylabel('Aceleración [g]')
title('Respuesta desp incial')
grid
nexttile
plot(t2,Respuesta_A_2)
xlabel('tiempo [s]')
ylabel('Aceleración [g]')
title('Respuesta impacto')
grid
nexttile
plot(t3,Respuesta_A_3)
xlabel('tiempo [s]')
ylabel('Aceleración [g]')
title('Respuesta amortiguamiento externo')
grid
nexttile
plot(t4,Respuesta_A_4)
xlabel('tiempo [s]')
ylabel('Aceleración [g]')
title('Respuesta resonancia')
grid

%Crear vector con l y p
l = [5.15 5.2 5.29 5.4 5.62 5.89 6.33]; %[cm]
p = [0 20 50 100 200 300 500]; %[gf]
coefficients = polyfit(l, p, 1); %Regresión lineal
k = coefficients(1)/10; %[kgf/m]
reg1 = coefficients(1).*p+l;

M = 0.635; %[kgF]
%Periodo teorico
w = sqrt(k/M);
T = (2*pi)/(w);

%Cortar el vector
Resp1_rec = Respuesta_A_1(2001:2801); %Desplazamiento inicial
Resp2_rec = Respuesta_A_2(5201:6001); %Velocidad inicial
Resp3_rec = Respuesta_A_3(101:461); %Amortiguamiento inicial
Resp4_rec = Respuesta_A_4(10001:10801); %Disipación forzante
Resp5_rec = Respuesta_A_4(2001:2801); %Resonacia

%Calculo de periodo

Fs = 200;

TF_g11 = fft(Resp1_rec);
L1 = length(Resp1_rec);
P2_g11 = abs(TF_g11/L1);
P1_g11 = P2_g11(1:L1/2+1);
P1_g11(2:end-1) = 2*P1_g11(2:end-1);
f_g11 = Fs*(0:(L1/2))/L1;
figure()
plot(f_g11,P1_g11)
grid
grid minor
xlabel('frecuencia [Hz]')
ylabel('Amplitud [g]')
title('Tranformada de Fourier Desplazamiento incial')
Max_P1_g11 = max(P1_g11);
indice_Max_P1_g11 = find(P1_g11 == Max_P1_g11);
f1 = f_g11(indice_Max_P1_g11);

TF_g12 = fft(Resp2_rec);
L2 = length(Resp2_rec);
P2_g12 = abs(TF_g12/L2);
P1_g12 = P2_g12(1:L2/2+1);
P1_g12(2:end-1) = 2*P1_g12(2:end-1);
f_g12 = Fs*(0:(L2/2))/L2;
figure()
plot(f_g12,P1_g12)
grid
grid minor
xlabel('frecuencia [Hz]')
ylabel('Amplitud [g]')
title('Tranformada de Fourier Impacto')
Max_P1_g12 = max(P1_g12);
indice_Max_P1_g12 = find(P1_g12 == Max_P1_g12);
f2 = f_g12(indice_Max_P1_g12);

TF_g13 = fft(Resp3_rec);
L3 = length(Resp3_rec);
P2_g13 = abs(TF_g13/L2);
P1_g13 = P2_g13(1:L3/2+1);
P1_g13(2:end-1) = 2*P1_g13(2:end-1);
f_g13 = Fs*(0:(L3/2))/L3;
figure()
plot(f_g13,P1_g13)
grid
grid minor
xlabel('frecuencia [Hz]')
ylabel('Amplitud [g]')
title('Tranformada de Fourier Amortiguamiento adicional')
Max_P1_g13 = max(P1_g13);
indice_Max_P1_g13 = find(P1_g13 == Max_P1_g13);
f3 = f_g13(indice_Max_P1_g13);

TF_g14 = fft(Resp4_rec);
L4 = length(Resp4_rec);
P2_g14 = abs(TF_g14/L4);
P1_g14 = P2_g14(1:L4/2+1);
P1_g14(2:end-1) = 2*P1_g14(2:end-1);
f_g14 = Fs*(0:(L4/2))/L4;
figure()
plot(f_g14,P1_g14)
grid
grid minor
xlabel('frecuencia [Hz]')
ylabel('Amplitud [g]')
title('Tranformada de Disipación Forzante')
Max_P1_g14 = max(P1_g14);
indice_Max_P1_g14 = find(P1_g14 == Max_P1_g14);
f4 = f_g14(indice_Max_P1_g14);

TF_g15 = fft(Resp5_rec);
L5 = length(Resp5_rec);
P2_g15 = abs(TF_g15/L5);
P1_g15 = P2_g15(1:L5/2+1);
P1_g15(2:end-1) = 2*P1_g15(2:end-1);
f_g15 = Fs*(0:(L5/2))/L5;
figure()
plot(f_g15,P1_g15)
grid
grid minor
xlabel('frecuencia [Hz]')
ylabel('Amplitud [g]')
title('Tranformada de Resonancia')
Max_P1_g15 = max(P1_g15);
indice_Max_P1_g15 = find(P1_g15 == Max_P1_g15);
f5 = f_g15(indice_Max_P1_g15);

%Calculo de beta
A_1_abs = abs(Resp1_rec);   
peaks_A1 = findpeaks(A_1_abs); %Peaks
num_maximos1 = linspace(1,length(peaks_A1),length(peaks_A1)); %Número de máximos

A_2_abs = abs(Resp2_rec);
peaks_A2 = findpeaks(A_2_abs);
num_maximos2 = linspace(1,length(peaks_A2),length(peaks_A2));

A_3_abs = abs(Resp3_rec);
peaks_A3 = findpeaks(A_3_abs);
num_maximos3 = linspace(1,length(peaks_A3),length(peaks_A3));

A_4_abs = abs(Resp4_rec);
peaks_A4 = findpeaks(A_4_abs);
num_maximos4 = linspace(1,length(peaks_A4),length(peaks_A4));

regresion1 = polyfit(num_maximos1, log(peaks_A1) , 1); 
beta1 = abs(regresion1(1))/pi;

regresion2 = polyfit(num_maximos2, log(peaks_A2) , 1); 
beta2 = abs(regresion2(1))/pi;

regresion3 = polyfit(num_maximos3, log(peaks_A3) , 1); 
beta3 = abs(regresion3(1))/pi;

regresion4 = polyfit(num_maximos4, log(peaks_A4) , 1); 
beta4 = abs(regresion4(1))/pi;

reg1 = regresion1(1).*num_maximos1+regresion1(2);
reg2 = regresion2(1).*num_maximos2+regresion2(2);
reg3 = regresion3(1).*num_maximos3+regresion3(2);
reg4 = regresion4(1).*num_maximos4+regresion4(2);

figure()
hold on
plot(num_maximos1,log(peaks_A1),'.')
plot(num_maximos1,reg1,'-')

figure()
hold on
plot(num_maximos2,log(peaks_A2),'.')
plot(num_maximos2,reg2,'-')

figure()
hold on
plot(num_maximos3,log(peaks_A3),'.')
plot(num_maximos3,reg3,'-')

figure()
hold on
plot(num_maximos4,log(peaks_A4),'.')
plot(num_maximos4,reg4,'-')

%Aceleración forzante vs aceleración de la respuesta
Resp_forz_rec = Respuesta_A_4(2001:2801);
figure
plot(Resp1_rec,Resp_forz_rec)
xlabel('Aceleración de respuesta')
ylabel('Aceleración del forzante')
grid

vo = peaks_A2(1)/(mean([f1,f2,f3,f4])*2*pi); %Velocidad inicial