%% GRÁFICAS CORRIENTE UV

clc
clear
close all
run Configuracion

%INPUTS

Tw=Tw; % Twindow (segundos)
Ts=Ts;  % Tstep  (segundos)


DirDataVelocidades= ['./Outputs/3-Velocidades/' NameV 'D_OutputVelocidades_Tw' num2str(Tw) '_Ts' num2str(Ts) '.mat']; %Datos Grillas Corrientes
load(DirDataVelocidades)

LC=load(['./Outputs/3-Velocidades/' NameV 'C_Output_Vydx' num2str(dx) 'dy' num2str(dy) '_Tw' num2str(Tw) 'Ts' num2str(Ts) '.mat']);

% Archivo de contiene resultados de algoritmo Chickadel
DirOutChick=['./Outputs/2-OutputChickadel/' NameV '_OutputChickadel_Tw' num2str(Tw) '_Ts' num2str(Ts) '.mat'];
load(DirOutChick)
%Tiempos
TCk=videoCurrentOut(200).InstVbar.t*24*3600; %Vector de tiempo Chikadel tiempo total OBS: Para todos las ubicaciones es el mismo vector
TCk=TCk-TCk(1);                            %tiempo inicializado en primer frame

dTstep=mean(diff(TCk)); %Ventana de tiempo, donde se calcula la velocidad media
TCk2=TCk+dTstep;        %Vector que contiene las ventanas de tiempos acumuladas
dNT=dTstep/0.5;


%Coordenadas y vector de tiempo asociado a cada frame rectificado 
DirDataF= ['./Inputs/DataFrames_' NameV];
LoadData=load(DirDataF);
DataFrames=LoadData.DataFrames;


% Ubicación imágenes rectificadas 
% (utilizado sólo para generar gráficos finales)
imageDirectory=['C:\Memoria2020\Resultados\Rectificados\' NameV '\MAT'];

% Lista con nombre y dirección de cada impagen
L=string(ls(imageDirectory));
L=L(3:end);
Nframes=length(L);    
for i=1:Nframes
    L(i)=strcat('FrameRec_',num2str(i),'.mat');
end


% Dimensiones Grilla
dim=size(x_grilla1);
Ni=dim(2);
Nj=dim(1);

%%
% Valores Ux truncados

MaxUx=3; %Valor absoluto

for kk=1:NT
    Ux_Trunc=GrillaUx(kk).Uxt2;
    for i=1:Ni
       for j=1:Nj
           if abs(Ux_Trunc(j,i)) > MaxUx     
               Ux_Trunc(j,i)=sign(Ux_Trunc(j,i))*MaxUx;
           end
       end
       end
    GrillaUx(kk).Ux_Trunc=Ux_Trunc;
end


%% Graficas Promedio

%Promedio y desviación estandar

for i=1:Ni
    for j=1:Nj
        for k =1:NT
            ListaUij(k)=GrillaUx(k).Uxt2(j,i);
            %ListaUij_Trunc(k)=GrillaUx(k).Ux_Trunc(j,i);
            %ListaVij(k)=GrillaVy(k).MV(j,i);
            ListaMagnitudij(k)=GrillaUx(k).Magnitud(j,i);
            
        end
        
        GCorrientes.UxProm2_mean(j,i)=mean(ListaUij,'omitnan');
        GCorrientes.UxProm2_st(j,i)=std(ListaUij,'omitnan');
        
        %GCorrientes.Vy_mean(j,i)=mean(ListaVij);
        %GCorrientes.Vy_st(j,i)=std(ListaVij);
        
        GCorrientes.MagnitudProm2_mean(j,i)=mean(ListaMagnitudij,'omitnan');
        GCorrientes.MagnitudProm2_st(j,i)=std(ListaMagnitudij,'omitnan');
        
        %GCorrientes.UxTrunc_mean(j,i)=mean(ListaUij_Trunc);
        %GCorrientes.UxTrunc_st(j,i)=std(ListaUij_Trunc);
        
    end
end



%% Grafica Vy

figura1=figure(1)
figura1.Position = [0 0 1400 1920]
U=zeros(length(y_g1),length(x_g1));

%GRÁFICA MATRIZ DE VELOCIDAD 
%---------------------------

subplot(1,2,1)

imagesc(LC.COutVy.x_list,LC.COutVy.y_list,GCorrientes.Vy_mean)
a = colorbar;
caxis([-1 1])
axis xy
xlim([80 300])
ylim([-270 260])
title('Mean Vy ')
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Velocidad [m/s]';
daspect([1 1 1])

hold on;



subplot(1,2,2)

imagesc(LC.COutVy.x_list,LC.COutVy.y_list,GCorrientes.Vy_st)
a = colorbar;
caxis([0 1.5])
axis xy
xlim([80 300])
ylim([-260 260])
title('Std Vy ')
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Velocidad [m/s]';
daspect([1 1 1])

sgtitle(['Mean Vy , Tw = ' num2str(Tw) ' s , Ts = ' num2str(Ts) ' s'])

saveas(figura1,['./Outputs/4-Resultados-Graficos/' NameV 'Vy_Tw' num2str(Tw) '_Ts' num2str(Ts) '.png'])


%% Grafica Ux (Promediando Vy)

figura2=figure(2)
U=zeros(length(y_g1),length(x_g1));

%GRÁFICA MATRIZ DE VELOCIDAD 
%---------------------------

subplot(1,1,1)

imagesc(x_g1,y_g1,GCorrientes.Ux2_mean)
a = colorbar;
caxis([-3 3])
axis xy
xlim([70 300])
ylim([-270 260])
title('Ux')
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Velocidad [m/s]';
daspect([1 1 1])

hold on;

sgtitle(['Ux (Mean Vy) , Tw = ' num2str(Tw) ' s , Ts = ' num2str(Ts) ' s'])
% 
% subplot(1,2,2)
% 
% imagesc(x_g1,y_g1,GCorrientes.Ux_st)
% a = colorbar;
% caxis([0 2])
% axis xy
% xlim([115 370])
% ylim([-270 250])
% title('Std Ux , Tw=256')
% xlabel('Cross-shore position [m]')
% ylabel('Alongshore position[m]')
% a.Label.String = 'Velocidad [m/s]';
% daspect([1 1 1])

figura2.Position = [0 0 1400 1920]
saveas(figura2,['./Outputs/4-Resultados-Graficos/' NameV '_PlotUx(MeanVy)_Tw' num2str(Tw) '_Ts' num2str(Ts) '.png'])


%% Grafica Ux (Promediando Ux)

figura3=figure(3)
U=zeros(length(y_g1),length(x_g1));

%GRÁFICA MATRIZ DE VELOCIDAD 
%---------------------------

subplot(1,2,1)

imagesc(x_g1,y_g1,GCorrientes.UxProm2_mean)
a = colorbar;
caxis([-3 3])
axis xy
xlim([75 300])
ylim([-270 260])
title('Mean Ux')
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Velocidad [m/s]';
daspect([1 1 1])

hold on;



subplot(1,2,2)

imagesc(x_g1,y_g1,GCorrientes.UxProm2_st)
a = colorbar;
caxis([0 2])
axis xy
xlim([75 300])
ylim([-270 260])
title('Std Ux' )
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Velocidad [m/s]';
daspect([1 1 1])

sgtitle(['Mean Ux , Tw = ' num2str(Tw) ' s , Ts = ' num2str(Ts) ' s'])
figura3.Position = [0 0 1400 1920]
saveas(figura3,['./Outputs/4-Resultados-Graficos/' NameV '_PlotUx(MeanUx)_Tw' num2str(Tw) '_Ts' num2str(Ts) '.png'])



%% Grafica Corrientes UV (PromediandoVy)

figura4=figure(4)



%GRÁFICA MATRIZ DE VELOCIDAD 
%---------------------------
subplot(1,1,1)
imagesc(x_g1,y_g1,GCorrientes.Magnitud_Mean)
a = colorbar;
caxis([-1 1])
axis xy
xlim([75 300])
ylim([-270 260])
title('Mean UV')
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Magnitud [m/s]';
daspect([1 1 1])

%Grafica Velocidad como Campo Vectorial
hold on;
quiver(x_grilla1,y_grilla1,GCorrientes.Ux2_mean,GCorrientes.VyF_mean,'r','LineWidth',1,'AutoScaleFactor',1.5)



% subplot(1,2,2)
% 
% imagesc(x_g1,y_g1,GCorrientes.Magnitud_st)
% a = colorbar;
% caxis([0 2])
% axis xy
% xlim([115 370])
% ylim([-270 250])
% title('Std UV ')
% xlabel('Cross-shore position [m]')
% ylabel('Alongshore position[m]')
% a.Label.String = 'Magnitud [m/s]';
% daspect([1 1 1])

sgtitle(['UV (MeanVy) , Tw = ' num2str(Tw) ' s , Ts = ' num2str(Ts) ' s'])
figura4.Position = [0 0 1400 1920]
saveas(figura4,['./Outputs/4-Resultados-Graficos/' NameV '_PlotUV(MeanVy)_Tw' num2str(Tw) '_Ts' num2str(Ts) '.png'])


%%
figure6=figure(6)
im6=imagesc(x_g1,y_g1,GCorrientes.Magnitud_Mean)
im6.AlphaData =0.7;
a = colorbar;
caxis([-1 1])
axis xy
xlim([65 310])
ylim([-270 260])
title('Mean UV')
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Magnitud [m/s]';
daspect([1 1 1])



f6=streamslice(x_grilla1,y_grilla1,GCorrientes.Ux2_mean,GCorrientes.VyF_mean,3,'r')

axis tight
set( f6, 'Color', 'k')
figure6.Position = [0 0 1400 1920]
saveas(figure6,['./Outputs/4-Resultados-Graficos/' NameV '_PlotUV(MeanVy-stream)_Tw' num2str(Tw) '_Ts' num2str(Ts) '.png'])
%% Grafica Corrientes UV (Promediando Ux)

figura5=figure(5)



%GRÁFICA MATRIZ DE VELOCIDAD 
%---------------------------
subplot(1,2,1)
imagesc(x_g1,y_g1,GCorrientes.MagnitudProm2_mean)
a = colorbar;
caxis([-1 1])
axis xy
xlim([65 310])
ylim([-270 260])
title('Mean UV')
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Magnitud [m/s]';
daspect([1 1 1])

%Grafica Velocidad como Campo Vectorial
hold on;
quiver(x_grilla1,y_grilla1,GCorrientes.UxProm2_mean,GCorrientes.VyF_mean,'r','LineWidth',1,'AutoScaleFactor',1.5)



subplot(1,2,2)

imagesc(x_g1,y_g1,GCorrientes.MagnitudProm2_st)
a = colorbar;
caxis([0 2])
axis xy
xlim([65 310])
ylim([-270 260])
title('Std UV ')
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Magnitud [m/s]';
daspect([1 1 1])

sgtitle(['UV (MeanUx) , Tw = ' num2str(Tw) ' s , Ts = ' num2str(Ts) ' s'])
figura5.Position = [0 0 1400 1920]
saveas(figura5,['./Outputs/4-Resultados-Graficos/' NameV '_PlotUV(MeanUx)_Tw' num2str(Tw) '_Ts' num2str(Ts) '.png'])



