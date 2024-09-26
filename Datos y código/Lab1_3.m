clc;close all;clear 

%Cargar archivos
g21 = load('g21.txt'); %Desplazamiento incial conocido
g22 = load('g22.txt'); %Golpe
g23 = load('g23.txt'); %Estructura con amortiguamiento
g24 = load('g24.txt'); %Resonancia 

%Constantes sensores
c_estructura = 1.029; 
c_forzante = 1.028;

%Vectores de tiempo de cada registro
t1 = g21(:,1);
t2 = g22(:,1);
t3 = g23(:,1);
t4 = g24(:,1);

%Respuesta aceleración sistema [g]
Respuesta_1 = g21(:,2)./c_estructura;
Respuesta_2 = g22(:,2)./c_estructura;
Respuesta_3 = g23(:,2)./c_estructura;
Respuesta_4 = g24(:,2)./c_estructura;

%Datos de la forzante [g]
Respuesta_for_1 = g21(:,3)./c_forzante;
Respuesta_for_2 = g22(:,3)./c_forzante;
Respuesta_for_3 = g23(:,3)./c_forzante;
Respuesta_for_4 = g24(:,3)./c_forzante;

%Gráfico respuestas aceleración
figure()
nexttile
plot(t1,Respuesta_1)
xlabel('Tiempo [s]')
ylabel('Aceleración [g]')
title('Respuesta con desplazamiento incial')
grid
nexttile
plot(t2,Respuesta_2)
xlabel('Tiempo [s]')
ylabel('Aceleración [g]')
title('Respuesta con golpe')
grid
nexttile
plot(t3,Respuesta_3)
xlabel('Tiempo [s]')
ylabel('Aceleración [g]')
title('Respuesta con amortiguamiento')
grid
nexttile
plot(t4,Respuesta_4)
xlabel('Tiempo [s]')
ylabel('Aceleración [g]')
title('Respuesta en resonancia')
grid

%Crear vector con l y p
l = [4.15 4.6 5.5 6.45 7.2 7.3 8.35]; %[cm]
p = [0 50 150 250 350 370 500]; %[gf]
coefficients = polyfit(l, p, 1); %Regresión lineal
k = coefficients(1)/10; %[kgf/m]
reg1 = coefficients(1).*p+l;

M = 0.635; %[kgF]
%Periodo teorico
w = sqrt(k/M);
T = (2*pi)/(w);

%Recortamos el vector de datos para calcular frecuencias naturales
Resp1_rec = Respuesta_1(2001:6001); %Desplazamiento inicial conocido
Resp2_rec = Respuesta_2(2001:6001); %Golpe
Resp3_rec = Respuesta_3(2543:3239); %Amortiguamiento
Resp4_rec = Respuesta_4(10001:20001); %Disipación forzante
Resp5_rec = Respuesta_4(2001:8001); %Resonacia

%Calculo de periodo
Fs = 200;

%Transformadas de fourier de las respuestas de la estructura
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
xlim([0 10])
ylabel('Amplitud [g]')
title('Tranformada de respuesta a un desplazamiento incial conocido')
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
xlim([0 7])
ylabel('Amplitud [g]')
title('Tranformada de respusta del golpe')
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
xlim([0 7])
ylabel('Amplitud [g]')
title('Tranformada con amortiguamiento adicional')
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
xlim([0 7])
ylabel('Amplitud [g]')
title('Tranformada de disipación forzante')
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
xlim([0 7])
ylabel('Amplitud [g]')
title('Tranformada de respuesta en resonancia')
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
Resp_forz_rec = Respuesta_4(2001:6001);
figure
plot(Resp1_rec,Resp_forz_rec)
xlabel('Aceleración de respuesta [g]')
ylabel('Aceleración del forzante [g]')
grid

vo = peaks_A2(1)/(mean([f1,f2,f3,f4])*2*pi); %Velocidad inicial