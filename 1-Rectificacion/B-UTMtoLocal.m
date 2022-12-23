% SISTEMA DE REFERENCIA LOCAL
clear
clc
% Ingresar nombre y coordenadas de los puntos de referencia [E N] en UTM
GCP_UTM.Name(1,:)="P1";
GCP_UTM.Name(2,:)="P2";
GCP_UTM.Name(3,:)="P3";
GCP_UTM.CoordUTM(1,:)=[256888.121 6289629.746]';
GCP_UTM.CoordUTM(2,:)=[256867.711 6289684.774]';
GCP_UTM.CoordUTM(3,:)=[256904.663 6289669.813]';


Ng=length(GCP_UTM.CoordUTM(:,1));

%%
%Ingresar punto de referencia [E N]
Re=[256884.448 6289741.483]';

%Ingresar el angulo de giro en grados (counterclockwise)
theta =150;

%Se genera la matriz de transformación
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];

%Se calculan los puntos referidos y se generan nuevas coordenadas´
for i=1:Ng;
    PuntosRef(i,:)= GCP_UTM.CoordUTM(i,:)-Re';
    RotP(i,:)=(R*PuntosRef(i,:)')'; %Nuevas Coordenadas
end
%%
%Se crea el archivo refPOINT

for i= 1:Ng;
    gcp(i).num=i;
    gcp(i).name=GCP_UTM.Name(i,:);
    gcp(i).x=RotP(i,1); 
    gcp(i).y=RotP(i,2);  
    gcp(i).z=0; %No debería afectar 
end  

clearvars -except gcp
save('./Outputs/refPOINT','gcp')
