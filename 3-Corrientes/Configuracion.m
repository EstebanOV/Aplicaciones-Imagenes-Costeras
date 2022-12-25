
%Vuelo a analizar
NameV = '15V3';

%Parámetros
dy = 10 ;  %[m] Longitud de los instrumentos vBar
dx = 10 ;  %[m] Separación instrumentos en dirección x
Tw = 32 ;  %[s] Duración ventana de análisis
Ts = 16 ;  %[s] Tiempo de muestreo
km = 0.25; %[1/m] Número de onda mínimo
P  = 0.9 ; %Límite para el nivel de confianza de ajuste al modelo
Lci= 0.2 ; %[m/s] Umbral para el intervalo de confianza de la velocidad media
Ir = 40  ; %Umbral para el rango de intensidad

%Límites grilla
Xlim=[70,370];
Ylim=[-270,260];