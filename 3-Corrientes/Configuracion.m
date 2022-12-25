
%Vuelo a analizar
NameV = '15V3';

%Par�metros
dy = 10 ;  %[m] Longitud de los instrumentos vBar
dx = 10 ;  %[m] Separaci�n instrumentos en direcci�n x
Tw = 32 ;  %[s] Duraci�n ventana de an�lisis
Ts = 16 ;  %[s] Tiempo de muestreo
km = 0.25; %[1/m] N�mero de onda m�nimo
P  = 0.9 ; %L�mite para el nivel de confianza de ajuste al modelo
Lci= 0.2 ; %[m/s] Umbral para el intervalo de confianza de la velocidad media
Ir = 40  ; %Umbral para el rango de intensidad

%L�mites grilla
Xlim=[70,370];
Ylim=[-270,260];