
clc
clear
close all
run Configuracion

%INPUTS

Tw=Tw; % Twindow (segundos)
Ts=Ts;  % Tstep  (segundos)

%Cargar resultados obtención Vy
DirVy=['./Outputs/3-Velocidades/' NameV 'C_Output_Vydx' num2str(dx) 'dy' num2str(dy) '_Tw' num2str(Tw) 'Ts' num2str(Ts) '.mat'];
load(DirVy)

x_g1=COutVy.x_list;     % Lista de X de la grilla Vy
y_g1=COutVy.y_list;     % Lista de Y de la grilla Vx y Vy?
GrillaVy=COutVy.GrillaVy;

%Cargar Batimetría
load('./Inputs/Bathy-5m-Las-Cruces.mat')

%% Compatibilización de Grillas
% (Deben estar en los mismos límites)

Xmax = max(bathy.xm); %295
Xmin = min(bathy.xm); %85
%Ymax = max(bathy.ym); Xmin = min(bathy.ym); %(Puede ser necesario)

XG = COutVy.GrillaX(:); YG = COutVy.GrillaY(:);

X = XG(XG>=Xmin & XG<=Xmax);
Y = YG(XG>=Xmin & XG<=Xmax);
x_g1=x_g1(x_g1>=Xmin) ; x_g1=x_g1(x_g1<=Xmax);

VM = GCorrientes.VyF_mean(:); VM = VM(XG>=Xmin & XG<=Xmax);
Vmean = reshape(VM,length(y_g1),length(x_g1)); 
GCorrientes.VyF_mean=Vmean;





for i=1:length(GrillaVy)
    GVi=GrillaVy(i).MV(:);
    GVi = GVi(XG>=Xmin & XG<=Xmax);
    GrillaVy2(i).MV=reshape(GVi,length(y_g1),length(x_g1)); 
end

GrillaVy=GrillaVy2;

% Generación Grilla para cálculo Vx (Grilla 2)
dx=mean(diff(x_g1));
dy=mean(diff(y_g1));
N=length(x_g1);
M=length(y_g1);

x_g2=(x_g1(1,1)-dx/2):dx:(x_g1(1,N)+dx/2);
y_g2=y_g1;


% Generación Grilla Batimetría
[x_grilla2,y_grilla2]=meshgrid(x_g2,y_g2);
Depth_g2=griddata(bathy.X(:),bathy.Y(:),bathy.h(:),x_grilla2,y_grilla2,'linear'); %Interpolación Lineal
Depth_g2 = inpaint_nans(Depth_g2);

%Grafica Batimetria
figure(1)

%Grafica profundidad cBathy (con filtro de Kalman)
imagesc(x_g1,y_g1,-Depth_g2); grid on
colormap parula
a = colorbar; a.Label.String = 'Profundidad [m]';
caxis([-10 0])
axis xy
xlim([80 300]); ylim([-265 255])
axis xy; xlabel('x (m)'); ylabel('y (m)');
title(['Profundidad cBathy'])
daspect([1 1 1])
hold on
contour(x_grilla2,y_grilla2,Depth_g2,'ShowText','on')

%% Cálculo Ux

NT=length(GrillaVy);
N2=length(x_g2);
M2=length(y_g2);
h=Depth_g2*(-1);       % Profundidad [m]
dx;
dy;


for kk=1:NT
kk;
V=GrillaVy(kk).MV;      % Velocidad Vy [m/s]
U=nan(M2,N2);          % Inicialización Matriz U 
U(:,1)=0;              % Condicion de Borde


for i=1:(N2-1)
   for j=1:(M2-1)
       
       hy1=mean([h(j,i) h(j,i+1)]);
       hy2=mean([h(j+1,i) h(j+1,i+1)]);
       hx1=mean([h(j,i) h(j+1,i)]);
       hx2=mean([h(j,i+1) h(j+1,i+1)]);
       V1=V(j,i);
       V2=V(j+1,i);
       U1=U(j,i);
       
       if ~isnan(U1)
           
           if (~isnan(V1)) && (~isnan(V2))
               
               U2=(1/hx2)*((dx/dy)*(hy1*V1-hy2*V2)+hx1*U1);
               U(j,i+1)=U2;
           end
           
           if (isnan(V1)) && (~isnan(V2))
               
               V1=0;
               U2=(1/hx2)*((dx/dy)*(hy1*V1-hy2*V2)+hx1*U1);
               U(j,i+1)=U2;
           end
           
           if (~isnan(V1)) && (isnan(V2))
               V2=0;
               U2=(1/hx2)*((dx/dy)*(hy1*V1-hy2*V2)+hx1*U1);
               U(j,i+1)=U2;
               
           end
           
           if (isnan(V1)) && (isnan(V2))
               
               U2=NaN;
               U(j,i+1)=U2;
           end
       end
       
       
       if isnan(U1)
           
           U1=0;   
           
           if (~isnan(V1)) && (~isnan(V2))
               
               U2=(1/hx2)*((dx/dy)*(hy1*V1-hy2*V2)+hx1*U1);
               U(j,i+1)=U2;
           end
           
           if (isnan(V1)) && (~isnan(V2))
               
               V1=0;
               U2=(1/hx2)*((dx/dy)*(hy1*V1-hy2*V2)+hx1*U1);
               U(j,i+1)=U2;
           end
           
           if (~isnan(V1)) && (isnan(V2))
               V2=0;
               U2=(1/hx2)*((dx/dy)*(hy1*V1-hy2*V2)+hx1*U1);
               U(j,i+1)=U2;
               
           end
           
           if (isnan(V1)) && (isnan(V2))
               
               U2=NaN;
               U(j,i+1)=U2;
           end
       end
   end
end







GrillaUx(kk).Uxt=U;
Uxt2=U;
Uxt2(:,1)=[];
GrillaUx(kk).Uxt2=Uxt2;


end



%% Calculo Ux promedio a partir de Vy promedio



Vmean;                 % Velocidad Vy [m/s]
U=nan(M2,N2);          % Inicialización Matriz U 
U(:,1)=0;              % Condicion de Borde


for i=1:(N2-1)
   for j=1:(M2-1)
       
       hy1=mean([h(j,i) h(j,i+1)]);
       hy2=mean([h(j+1,i) h(j+1,i+1)]);
       hx1=mean([h(j,i) h(j+1,i)]);
       hx2=mean([h(j,i+1) h(j+1,i+1)]);
       V1=Vmean(j,i);
       V2=Vmean(j+1,i);
       U1=U(j,i);
       
       if ~isnan(U1)
           
           if (~isnan(V1)) && (~isnan(V2))
               
               U2=(1/hx2)*((dx/dy)*(hy1*V1-hy2*V2)+hx1*U1);
               U(j,i+1)=U2;
           end
           
           if (isnan(V1)) && (~isnan(V2))
               
               V1=0;
               U2=(1/hx2)*((dx/dy)*(hy1*V1-hy2*V2)+hx1*U1);
               U(j,i+1)=U2;
           end
           
           if (~isnan(V1)) && (isnan(V2))
               V2=0;
               U2=(1/hx2)*((dx/dy)*(hy1*V1-hy2*V2)+hx1*U1);
               U(j,i+1)=U2;
               
           end
           
           if (isnan(V1)) && (isnan(V2))
               
               U2=NaN;
               U(j,i+1)=U2;
           end
       end
       
       
       if isnan(U1)
           
           U1=0;   
           
           if (~isnan(V1)) && (~isnan(V2))
               
               U2=(1/hx2)*((dx/dy)*(hy1*V1-hy2*V2)+hx1*U1);
               U(j,i+1)=U2;
           end
           
           if (isnan(V1)) && (~isnan(V2))
               
               V1=0;
               U2=(1/hx2)*((dx/dy)*(hy1*V1-hy2*V2)+hx1*U1);
               U(j,i+1)=U2;
           end
           
           if (~isnan(V1)) && (isnan(V2))
               V2=0;
               U2=(1/hx2)*((dx/dy)*(hy1*V1-hy2*V2)+hx1*U1);
               U(j,i+1)=U2;
               
           end
           
           if (isnan(V1)) && (isnan(V2))
               
               U2=NaN;
               U(j,i+1)=U2;
           end
       end
   end
end








Uxt2=U;
Uxt2(:,1)=[];

GCorrientes.Ux_mean=U;
GCorrientes.Ux2_mean=Uxt2;


%% CONTRUCCIÓN FRAMES PROMEDIADOS SEGÚN VENTANA DE TIEMPO 
%----------------------------------------------------------------

% Archivo de contiene resultados de algoritmo Chickadel
DirOutChick=['./Outputs/2-OutputChickadel/' NameV '_OutputChickadel_Tw' num2str(Tw) '_Ts' num2str(Ts) '.mat'];
load(DirOutChick)
%Tiempos
TCk=videoCurrentOut(200).InstVbar.t*24*3600; %Vector de tiempo Chikadel tiempo total OBS: Para todos las ubicaciones es el mismo vector
TCk=TCk-TCk(1);                            %tiempo inicializado en primer frame

dTstep=mean(diff(TCk)); %Ventana de tiempo, donde se calcula la velocidad media
TCk2=TCk+dTstep;        %Vector que contiene las ventanas de tiempos acumuladas
dNT=dTstep/0.5;


% Ubicación imágenes rectificadas 
% (utilizado sólo para generar gráficos finales)
imageDirectory=['C:\Memoria2020\Resultados\Rectificados\' NameV '\MAT\']; %Carpeta donde se ubican los frames

% Lista con nombre y dirección de cada imagen
L=string(ls(imageDirectory));
L=L(3:end);
Nframes=length(L);    
for i=1:Nframes
    L(i)=strcat('FrameRec_',num2str(i),'.mat');
end

% Frames promediados quedan almacenados en Variable Iprom

c1=0;
for i=1:length(TCk)
    iTimex=double(load(strcat(imageDirectory, L(200))).Irect).*0;
    c1=int16(1+dNT*(i-1));
    c2=int16(dNT*i);
    for j=c1:c2
        iTimex=iTimex+double(load(strcat(imageDirectory, L(j))).Irect);
        if j==c2
        IProm{i}=uint8(iTimex./dNT);   %Cálculo promedio 
        end
    end
end


%Coordenadas y vector de tiempo asociado a cada frame rectificado 
DirDataF= ['./Inputs/DataFrames_' NameV];
LoadData=load(DirDataF);
DataFrames=LoadData.DataFrames;







%% Gráficas Finales

TCk=COutVy.TCk;
dTstep=mean(diff(TCk)); %Ventana de tiempo, donde se calcula la velocidad media
TCk2=TCk+dTstep;        %Vector que contiene las ventanas de tiempos acumuladas

%Grafica Velocidad Ux
% figura10=figure(10)
% 
% for i=1:NT
% %GRÁFICA MATRIZ DE VELOCIDAD 
% %---------------------------
% 
% subplot(1,2,2)
% imagesc(x_g1,y_g1,GrillaUx(i).Uxt2) %Es correcta la utilización de esas coordendas? x_g1 y x_g2?
% a = colorbar;
% caxis([-3 3])
% axis xy
% xlim([115 300])
% ylim([-270 260])
% title('Ux , T = '+string(TCk2(i))+' s')
% xlabel('Cross-shore position [m]')
% ylabel('Alongshore position[m]')
% a.Label.String = 'Velocidad Ux [m/s]';
% daspect([1 1 1])
% hold on;
% 
% % GRÁFICA FRAMES PROMEDIADOS
% %----------------------------
% 
% subplot(1,2,1)
% 
% image(DataFrames.Coord.x(1,:),DataFrames.Coord.y(1,:),IProm{1,i});
% daspect([1 1 1])
% xlabel('Cross-shore position [m]')
% ylabel('Alongshore position[m]')
% axis xy
% xlim([115 300])
% ylim([-270 250])
% title('Imágenes Promediadas: \Deltat = 32 s')
% 
% %saveas(figura10,['./Resultados/Ux/Ux2_T' num2str(TCk2(i)) 's.png'])
% pause(0.01)
% 
% end
%%
% Valores Ux truncados

MaxUx=3; %Valor absoluto

for kk=1:NT
    Ux_Trunc=GrillaUx(kk).Uxt2;
    for i=1:(N2-1)
       for j=1:(M2)
           if abs(Ux_Trunc(j,i)) > MaxUx     
               Ux_Trunc(j,i)=sign(Ux_Trunc(j,i))*MaxUx;
           end
       end
       end
    GrillaUx(kk).Ux_Trunc=Ux_Trunc;
end

%% Magnitud

for kk=1:NT
    
    for i=1:(N2-1)
       for j=1:(M2)
           signoUx=sign(GrillaUx(kk).Ux_Trunc(j,i));
           signoVy=sign(GrillaVy(kk).MV(j,i));
           
           if signoVy<0
               
               Magnitud(j,i)=sqrt(GrillaUx(kk).Ux_Trunc(j,i)^2+GrillaVy(kk).MV(j,i)^2)*(-1);
               
           else
               
               Magnitud(j,i)=sqrt(GrillaUx(kk).Ux_Trunc(j,i)^2+GrillaVy(kk).MV(j,i)^2)*signoUx*signoVy;
               
           end
                 
       end
       end
    GrillaUx(kk).Magnitud=Magnitud;
end

%% Magnitud con Velocidades Promediadas


for i=1:(N2-1)
       for j=1:(M2)
           signoUx=sign(GCorrientes.Ux2_mean(j,i));
           signoVy=sign(GCorrientes.VyF_mean(j,i));
           
           if (signoUx<0) && (signoVy<0)
               
               MagnitudM(j,i)=sqrt(GCorrientes.Ux2_mean(j,i)^2+GCorrientes.VyF_mean(j,i)^2)*(-1);
               
           else
               
               MagnitudM(j,i)=sqrt(GCorrientes.Ux2_mean(j,i)^2+GCorrientes.VyF_mean(j,i)^2)*signoUx*signoVy;
               
           end
                 
       end
       end
   
GCorrientes.Magnitud_Mean=MagnitudM;




%%
% Grafica Corrientes
[x_grilla1,y_grilla1]=meshgrid(x_g1,y_g1);
% for i=1:NT
% figure(11)
% 
% %GRÁFICA MATRIZ DE VELOCIDAD 
% %---------------------------
% 
% 
% imagesc(x_g1,y_g1,GrillaUx(i).Magnitud)
% a = colorbar;
% caxis([-1 1])
% axis xy
% xlim([115 370])
% ylim([-270 250])
% %title('T = '+string(TCk2(i))+' s')
% xlabel('Cross-shore position [m]')
% ylabel('Alongshore position[m]')
% a.Label.String = 'Magnitud [m/s]';
% daspect([1 1 1])
% 
% 
% 
% hold on;
% quiver(x_grilla1,y_grilla1,GrillaUx(i).Ux_Trunc,GrillaVy(i).MV,'r','LineWidth',1,'AutoScaleFactor',1.5)
% pause(0.2)
% end


%% Sección 6: Output/Saving

% Output Name
oname=[NameV 'D_OutputVelocidades_Tw' num2str(Tw) '_Ts' num2str(Ts) '.mat'];

% OutPut Directory
odir='./Outputs/3-Velocidades';

% Save file
save([odir '/' oname ],'GrillaVy','GrillaUx','x_g1','y_g1','NT','x_grilla1','y_grilla1','GCorrientes')



