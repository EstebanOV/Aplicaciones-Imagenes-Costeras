%% OBTENCIÓN DE VELOCIDADES A PARTIR DE IMÁGENES DE DRON RECTIFICADAS

clc
clear
close all
addpath(genpath('.')) %Añadir carpetas

%% Parametros

edit Configuracion.m

%% Creacion instrumentos

tic
run A_CreacionInstrumento.m 
Tiempos.TA = toc;
save(['Set' NameV],'Tiempos')

% run Configuracion.m
% Tiempos.TA = 0;
% save(['Set' NameV],'Tiempos')


% Algoritmo de Chickadel
tic
run B_CorrienteChickadel.m

load(['Set' NameV])
Tiempos.TB = toc;
save(['Set' NameV],'Tiempos')

% Arreglo Grilla
tic
run C_ArregloGrillaV.m

load(['Set' NameV])
Tiempos.TC = toc;
save(['Set' NameV],'Tiempos')



% Obtención Ux
tic
run D_ObtencionUx.m

load(['Set' NameV])
Tiempos.TD = toc;
Tiempos.TT=(Tiempos.TA+Tiempos.TB+Tiempos.TC+Tiempos.TD)/60;
save(['Set' NameV],'Tiempos')


% Graficos Corrientes

run E_GraficaCorrientesUV.m









