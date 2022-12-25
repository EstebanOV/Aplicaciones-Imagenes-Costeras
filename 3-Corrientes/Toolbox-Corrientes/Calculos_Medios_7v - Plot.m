%% ARREGLO MATRICES DE VELOCIDADES
%==================================

clc
clear
close all
run Configuracion

%INPUTS

Tw=Tw; % Twindow (en puntos, no segundos)
Ts=Ts; % Tstep   (en puntos, no segundos)
                                              
% Archivos de contienen resultados Chickadel (Outputs Script B)
B12V1 = load(['./Outputs/2-OutputChickadel/12V1_OutputChickadel_Tw' num2str(Tw) '_Ts' num2str(Ts) '.mat']);
B13V1 = load(['./Outputs/2-OutputChickadel/13V1_OutputChickadel_Tw' num2str(Tw) '_Ts' num2str(Ts) '.mat']);
B13V2 = load(['./Outputs/2-OutputChickadel/13V2_OutputChickadel_Tw' num2str(Tw) '_Ts' num2str(Ts) '.mat']);
B13V3 = load(['./Outputs/2-OutputChickadel/13V3_OutputChickadel_Tw' num2str(Tw) '_Ts' num2str(Ts) '.mat']);
B15V1 = load(['./Outputs/2-OutputChickadel/15V1_OutputChickadel_Tw' num2str(Tw) '_Ts' num2str(Ts) '.mat']);
B15V2 = load(['./Outputs/2-OutputChickadel/15V2_OutputChickadel_Tw' num2str(Tw) '_Ts' num2str(Ts) '.mat']);
B15V3 = load(['./Outputs/2-OutputChickadel/15V3_OutputChickadel_Tw' num2str(Tw) '_Ts' num2str(Ts) '.mat']);

%Archivos arreglos velocidad Vy (Outputs Scripts C)
C12V1 = load(['./Outputs/3-Velocidades/12V1C_Output_Vydx' num2str(dx) 'dy' num2str(dy) '_Tw' num2str(Tw) 'Ts' num2str(Ts) '.mat']);
C13V1 = load(['./Outputs/3-Velocidades/13V1C_Output_Vydx' num2str(dx) 'dy' num2str(dy) '_Tw' num2str(Tw) 'Ts' num2str(Ts) '.mat']);
C13V2 = load(['./Outputs/3-Velocidades/13V2C_Output_Vydx' num2str(dx) 'dy' num2str(dy) '_Tw' num2str(Tw) 'Ts' num2str(Ts) '.mat']);
C13V3 = load(['./Outputs/3-Velocidades/13V3C_Output_Vydx' num2str(dx) 'dy' num2str(dy) '_Tw' num2str(Tw) 'Ts' num2str(Ts) '.mat']);
C15V1 = load(['./Outputs/3-Velocidades/15V1C_Output_Vydx' num2str(dx) 'dy' num2str(dy) '_Tw' num2str(Tw) 'Ts' num2str(Ts) '.mat']);
C15V2 = load(['./Outputs/3-Velocidades/15V2C_Output_Vydx' num2str(dx) 'dy' num2str(dy) '_Tw' num2str(Tw) 'Ts' num2str(Ts) '.mat']);
C15V3 = load(['./Outputs/3-Velocidades/15V3C_Output_Vydx' num2str(dx) 'dy' num2str(dy) '_Tw' num2str(Tw) 'Ts' num2str(Ts) '.mat']);

%Cargar Batimetría
load('./Inputs/Bathy-5m-Las-Cruces.mat')

%Velocidades
V1=C12V1.GCorrientes.Vy_mean;
V2=C13V1.GCorrientes.Vy_mean;
V3=C13V2.GCorrientes.Vy_mean;
V4=C13V3.GCorrientes.Vy_mean;
V5=C15V1.GCorrientes.Vy_mean;
V6=C15V2.GCorrientes.Vy_mean;
V7=C15V3.GCorrientes.Vy_mean;

%Coordenadas
x_g1 = C12V1.COutVy.x_list; 
y_g1 = C12V1.COutVy.y_list; 

%Grilla Vy
[x_grilla1,y_grilla1]=meshgrid(x_g1,y_g1);

% Información del instrumento
Name='12V1'; %Todos son iguales
DirIns=['./Outputs/1-Instrumentos/' Name '_Instrumento_Ly' num2str(dy) '_dx' num2str(dx) '.mat']; 
LoadIns=load(DirIns);Instr=LoadIns.instr;



%% %Filtro Velocidades
% V1(V1>1.5)=1.5; V1(V1<-1.5)=-1.5;  %Para otros datos puede ser útil

% Promedio Velocidades
VM=reshape(nanmean([V1(:)';V2(:)';V3(:)';V4(:)';V5(:)';V6(:)';V7(:)';]),size(V1));
VMstd=reshape(nanstd([V1(:)';V2(:)';V3(:)';V4(:)';V5(:)';V6(:)';V7(:)';]),size(V1));
MeanSTD=nanmean(VMstd(:));
figura1=figure(1);
subplot(1,2,1)
aap=imagesc(x_g1,y_g1,VM)
a = colorbar;
caxis([-1 1]);axis xy;daspect([1 1 1]);
xlim([80 300]);ylim([-270 260]);
title('Mean Vy')
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Velocidad [m/s]';
set(aap, 'AlphaData' , 1-isnan(VM))

subplot(1,2,2)
bbp=imagesc(x_g1,y_g1,VMstd)
a = colorbar;
caxis([0 1]);axis xy;daspect([1 1 1]);
xlim([80 300]);ylim([-270 260]);
title(['Std , Mean(Std) = ' num2str(round(MeanSTD,2)) ' m/s'])
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Velocidad [m/s]';
set(bbp, 'AlphaData' , 1-isnan(VMstd))

sgtitle(['Velocidad Media (12 - 13 - 15 nov)'])
figura1.Position = [0 0 1400 1920]
saveas(figura1,['./Outputs/5-Graficas-Promedio7v/VyMean-Std-7v.png'])
%%
%Filtro de Velocidades Altas
VMn = VM; VMn(abs(VMn)>0.6) = nan; %VMn(VMn>0.4) = nan;  0.56  Defina umbral si lo desea 

%Interpolación en datos nan
VMi= inpaint_nans(VMn);

%Filtro Gaussiano
%VMg = VMi ; %Cambie si quiere filtro
VMg=imgaussfilt(VMi); %filtra la matriz con un kernel de suavizado gaussiano 2D con una desviación estándar de 0,5

% Plot Velocidades
figura2=figure(2);
subplot(1,4,1)
aap=imagesc(x_g1,y_g1,VM)
a = colorbar;
caxis([-1 1]);axis xy;daspect([1 1 1]);
xlim([80 300]);ylim([-270 260]);
title('Mean Vy ')
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Velocidad [m/s]';
set(aap, 'AlphaData' , 1-isnan(VMn))

subplot(1,4,2)
bbp=imagesc(x_g1,y_g1,VMi)
a = colorbar;
caxis([-1 1]);axis xy;daspect([1 1 1]);
xlim([80 300]);ylim([-270 260]);
title('Filtro Umbral e Interpolación')
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Velocidad [m/s]';
set(bbp, 'AlphaData' , 1-isnan(VMi))

subplot(1,4,3)
ccp=imagesc(x_g1,y_g1,VMg)
a = colorbar;
caxis([-1 1]);axis xy;daspect([1 1 1]);
xlim([80 300]);ylim([-270 260]);
title('Filtro Gaussiano')
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Velocidad [m/s]';
% set(bbp, 'AlphaData' , 1-isnan(VMi))

subplot(1,4,4)
ddp=imagesc(x_g1,y_g1,abs(VMi-VMg))
a = colorbar;
caxis([0 0.5]);axis xy;daspect([1 1 1]);
xlim([80 300]);ylim([-270 260]);
title('Diff (Vint-Vfgauss) ')
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Velocidad [m/s]';
% set(bbp, 'AlphaData' , 1-isnan(VMi))

sgtitle(['Interpolación y Filtrado (12 - 13 - 15 nov)'])
figura2.Position = [0 0 1400 1920]
saveas(figura2,['./Outputs/5-Graficas-Promedio7v/Vy-Int_Filt-7v.png'])

%%  Cálculos previos a Ux

% Generación Grilla para cálculo Ux (Grilla 2)
dx = mean(diff(x_g1)); dy = mean(diff(y_g1));
N = length(x_g1); M = length(y_g1);

x_g2=(x_g1(1,1)-dx/2):dx:(x_g1(1,N)+dx/2);
y_g2=y_g1;

% Generación Grilla Batimetría
[x_grilla2,y_grilla2]=meshgrid(x_g2,y_g2);
Depth_g2=griddata(bathy.X(:),bathy.Y(:),bathy.h(:),x_grilla2,y_grilla2,'linear'); %Interpolación Lineal
Depth_g2 = inpaint_nans(Depth_g2);

% %Grafica Batimetria
% figure(3)
% 
% %Grafica profundidad cBathy (con filtro de Kalman)
% imagesc(x_g1,y_g1,-Depth_g2); grid on
% colormap parula
% a = colorbar; a.Label.String = 'Profundidad [m]';
% caxis([-10 0])
% axis xy
% xlim([80 300]); ylim([-265 255])
% axis xy; xlabel('x (m)'); ylabel('y (m)');
% title(['Profundidad cBathy'])
% daspect([1 1 1])
% hold on
% contour(x_grilla2,y_grilla2,Depth_g2,'ShowText','on')


%% CÁLCULOS Y GRAFICAS FINALES


%Vectores Ux y Vy
VM = VMi; UM = Ux; 

%Grilla Plot
xplot=x_g1; xplot=xplot+dx/2; xplot(length(xplot))=[];
yplot=y_g1; yplot=yplot+dy/2; yplot(length(yplot))=[];
[Xplot,Yplot]=meshgrid(xplot,yplot);

%Compatibilización de grillas
V_p = griddata(x_grilla1(:),y_grilla1(:),VM(:),Xplot,Yplot,'linear'); %Interpolación Lineal
U_p = griddata(xu_g(:),yu_g(:),UM(:),Xplot,Yplot,'linear'); %Interpolación Lineal

%Magnitud
M_p = sqrt(U_p(:)'.^2+V_p(:)'.^2)'; M_p = reshape(M_p,size(Xplot));

%Vectores Directores
DIR_p=[U_p(:)';V_p(:)']./sqrt(U_p(:)'.^2+V_p(:)'.^2);

%Filtro Velocidades > 1.5 m/s 
% Mag2 = Mag; Mag2(Mag2>1.5) = 1.5; 
% UM = DIR(1,:)'.*Mag2(:) ; VM = DIR(2,:)'.*Mag2(:) ; 
% UM = reshape(UM,size(VMi)) ; VM = reshape(VM,size(VMi));


%Grafica Vy - Ux
figura4=figure(4)
subplot(1,2,1)
imagesc(x_g1,y_g1,VM)
a = colorbar;caxis([-1 1]);axis xy;daspect([1 1 1])
xlim([75 300]);ylim([-270 260])
title('Vy')
xlabel('Cross-shore position [m]');ylabel('Alongshore position[m]')
a.Label.String = 'Magnitud [m/s]';

subplot(1,2,2)
imagesc(xu,yu,UM)
a = colorbar;caxis([-2 2]);axis xy;daspect([1 1 1])
xlim([75 300]);ylim([-270 260])
title('Ux')
xlabel('Cross-shore position [m]');ylabel('Alongshore position[m]')
a.Label.String = 'Magnitud [m/s]';
sgtitle(['Campo de Velocidades (12 - 13 - 15 nov)'])
figura4.Position = [0 0 1400 1920]
saveas(figura4,['./Outputs/5-Graficas-Promedio7v/Vy-Ux-7v.png'])

%Grafica Corrientes

MagPlot=M_p; %Mag.*sign(VM); %Magnitud
figura5=figure(5)
subplot(1,2,1)
imagesc(xplot,yplot,MagPlot)
a = colorbar;
caxis([0 1.2])
axis xy
xlim([75 300])
ylim([-270 260])
title('Corrientes UV')
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Magnitud [m/s]';
daspect([1 1 1])
%Grafica Velocidad como Campo Vectorial
hold on;
quiver(Xplot,Yplot,U_p,V_p,'r','LineWidth',1,'AutoScaleFactor',1.5)



subplot(1,2,2)
im7=imagesc(xplot,yplot,MagPlot)
im7.AlphaData =0.7;
a = colorbar;
caxis([0 1.2])
axis xy
xlim([75 300])
ylim([-270 260])
title('Corrientes UV')
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Magnitud [m/s]';
daspect([1 1 1])
f7=streamslice(Xplot,Yplot,U_p,V_p,3,'r');
axis tight ; set( f7, 'Color', 'k')

sgtitle(['Campo de Velocidades (12 - 13 - 15 nov)'])
figura5.Position = [0 0 1400 1920]
saveas(figura5,['./Outputs/5-Graficas-Promedio7v/Corrientes-7v.png'])

%% Output/Saving
%Variables
VM;UM;
xu;yu;
Xu = xu_g; Yu = yu_g;
xv = x_g1; yv = y_g1;
Xv = x_grilla1; Yv = y_grilla1;
VarPlot.xplot = xplot;
VarPlot.yplot = yplot;
VarPlot.Xplot = Xplot;
VarPlot.Yplot = Yplot;
VarPlot.U_p = U_p;
VarPlot.V_p = V_p;
VarPlot.M_p = M_p;
VarPlot.DIR_p = DIR_p;

% Output Name
Name='7vuelos';
oname=[Name '_CalculosMedios'];

% OutPut Directory
odir='./Outputs/6-Valores-Medios';

% Save file
save([odir '/' oname ],'VM','UM','xu','yu','Xu','Yu','xv','yv','Xv','Yv','VarPlot')
