%% OBTENCIÓN DE VELOCIDADES A PARTIR DE IMÁGENES DE DRON RECTIFICADAS

clc
clear
close all
addpath(genpath('.')) %Añadir carpetas
run Configuracion

%% Calculos Medios 7 vuelos

run F_Calculos_Medios_7v.m

%% Calculos Medios 15 Nov

run G_Calculos_Medios_15nov.m

%% Comparación Método Optico - Drifters

run H_Comparacion.m

close all