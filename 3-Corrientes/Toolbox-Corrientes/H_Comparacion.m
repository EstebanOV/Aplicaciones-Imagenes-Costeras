%% COMPARACION DATOS DRIFTERS CON CORRIENTES OBTENIDAS DE CHICKADEL

clc
clear
close all
run Configuracion

%INPUTS

%Velocidades Drifters
LDrift = load('./Inputs/Drifters/Velocidades_Medias_14-15Nov-Correccion.mat');

%Velocidades Chickadel - Continuidad
LChick = load('./Outputs/6-Valores-Medios/7vuelos_CalculosMedios.mat');

%Cargar Batimetría
load('./Inputs/Bathy-5m-Las-Cruces.mat')

%Definición celdas (Utilizadas en Drifters)
deltax=15; %[m]     tamaño celdas
deltay=15; %[m]     tamaño celdas
xg=[85:15:300];   % Límites grilla
yg=[-300:15:300]; % Límites grilla

%Grila X-Y
xplot=xg+deltax*0.5;xplot(length(xplot))=[];
yplot=yg+deltay*0.5;yplot(length(yplot))=[];
[XG,YG]=meshgrid(xplot,yplot);

%Coeficientes de ajuste
kx = -10;
ky = -10;

%Grilla U
UG = griddata(LChick.Xu(:)+kx,LChick.Yu(:)+ky,LChick.UM(:),XG,YG,'linear'); %Interpolación Lineal

%Grilla V
VG = griddata(LChick.Xv(:)+kx,LChick.Yv(:)+ky,LChick.VM(:),XG,YG,'linear'); %Interpolación Lineal

%Magnitud
MG = sqrt(UG(:)'.^2+VG(:)'.^2)'; MG = reshape(MG,size(XG));

%Direccicion
DG = [UG(:)';VG(:)']./sqrt((UG(:)').^2+(VG(:)').^2);
DIR1 = DG(1,:)'; DIR2 = DG(2,:)';
Dg1 = reshape(DIR1,size(XG)); 
Dg2 = reshape(DIR2,size(XG));
Dg1(isnan(LDrift.V_G)) = nan;
Dg2(isnan(LDrift.V_G)) = nan;


%% Grafica Magnitud Velocidad de Corriente

figura1 = figure(1);
clf

subplot(1,3,1)
h1=imagesc(xplot,yplot,MG);
colormap parula;a = colorbar; daspect([1 1 1])
caxis([0 1]);axis xy
xlim([100 200]);ylim([-100 150])
set(h1, 'AlphaData' , 1-isnan(LDrift.V_G))
hold on
quiver(XG(:),YG(:),Dg1(:),Dg2(:),0.3,'k')

title({'Chickadel - Continuidad'})
xlabel('x [m]');ylabel('y [m]')
a.Label.String = 'Velocidad [m/s]';
%Curvas de nivel
levels=[1 2 3 4];
[Cc,hc] = contour(bathy.X,bathy.Y,bathy.h,levels);
clabel(Cc,hc)

subplot(1,3,2)
h2=imagesc(xplot,yplot,LDrift.V_G);
colormap jet;b = colorbar; daspect([1 1 1])
caxis([0 1]);axis xy
xlim([100 200]);ylim([-100 150])
set(h2, 'AlphaData' , 1-isnan(LDrift.V_G))
hold on
quiver(LDrift.Xplot(:),LDrift.Yplot(:),LDrift.DIR1(:),LDrift.DIR2(:),0.3,'k')

title({'Drifters'})
xlabel('x [m]');ylabel('y [m]')
b.Label.String = 'Velocidad [m/s]';
%Curvas de nivel
levels=[1 2 3 4];
[Cc,hc] = contour(bathy.X,bathy.Y,bathy.h,levels);
clabel(Cc,hc)

Merr=nanmean(abs(MG(:)-LDrift.V_G(:)));
Mstd=nanstd(abs(MG(:)-LDrift.V_G(:)));
DiffM_p=abs(((MG-LDrift.V_G)./LDrift.V_G)*100); %Diferencia Porcentual
DiffM=abs(MG-LDrift.V_G); %Diferencia

subplot(1,3,3)
h3=imagesc(xplot,yplot,DiffM_p);
colormap jet;c = colorbar; daspect([1 1 1])
caxis([0 100]);axis xy
xlim([100 200]);ylim([-100 150])
set(h3, 'AlphaData' , 1-isnan(LDrift.V_G))
hold on
title(['Mean: ' num2str(round(Merr,3)) '0 Std: ' num2str(round(Mstd,3))])
xlabel('x [m]');ylabel('y [m]')
c.Label.String = 'Diferencia Porcentual [%]';

%Curvas de nivel
levels=[1 2 3 4];
[Cc,hc] = contour(bathy.X,bathy.Y,bathy.h,levels);
clabel(Cc,hc)

sgtitle(['Magnitud Velocidad'])

figura1.Position = [0 0 1400 1920]
saveas(figura1,['./Outputs/7-Comparacion/CompsUV.png'])

DiffMV=abs(MG(:)-LDrift.V_G(:));
Comp.Mag.diffMVmean=Merr;
Comp.Mag.diffMVstd=Mstd;
Comp.Mag.diffMVrms=sqrt(mean(DiffMV .* DiffMV, 'omitnan'));

%% Correlación entre mediciones in situ y mediciones ópticas
Mi = LDrift.V_G(:);                 % Medición in situ (Drifters) 
Ms = MG(:);     % Medición óptica (Chickadel-Continuidad)

Mi(isnan(Ms)) = nan;
Ms(isnan(Mi)) = nan;
Mi(isnan(Mi)) = [];
Ms(isnan(Ms)) = [];

%Regresión lineal: Ym = beta0 + beta1 * Xm

Xm = [ones(length(Mi),1) Mi];
b = Xm\Ms;
Ym = Xm*b;
R2 = 1 - sum((Ms - Ym).^2)/sum((Ms - mean(Ms)).^2);
Comp.Mag.R2=R2;

x_line=[-0.2:0.1:1.5]';
X_line = [ones(length(x_line),1) x_line];
Y_line = X_line*b;

figura4 = figure(4);
clf
scatter(Mi,Ms)
hold on
pv = plot(x_line,Y_line)
pv(1).LineWidth = 1.5;
xlim([-0.1 1.4]);ylim([-0.3 1.2])

xlabel('Mi [m/s] (mediciones drifters)')
ylabel('Ms [m/s] (método óptico)')

legend({'Data','Regresión Lineal'})
daspect([1 1 1])

txt = [{['Vs = ' num2str(round(b(2),3)) 'Vi + ' num2str(round(b(1),3))]},{['R^2 = ' num2str(round(R2,3))]}];
text(1,0.6,txt)

figura4.Position = [0 0 1400 1920]
saveas(figura4,['./Outputs/7-Comparacion/RegresionUV.png'])



%% Vy

figura2 = figure(2)
clf

subplot(1,3,1)
h1=imagesc(xplot,yplot,VG);
colormap jet;a = colorbar; daspect([1 1 1])
caxis([-1 1]);axis xy
xlim([100 200]);ylim([-100 150])
set(h1, 'AlphaData' , 1-isnan(LDrift.V_G))
title({'Chickadel - Continuidad'})
xlabel('x [m]');ylabel('y [m]')
a.Label.String = 'Velocidad [m/s]';
hold on
%Curvas de nivel
levels=[1 2 3 4];
[Cc,hc] = contour(bathy.X,bathy.Y,bathy.h,levels);
clabel(Cc,hc)


subplot(1,3,2)

%Vd=reshape(LDrift.V,flip(size(XG)))';
Vd=LDrift.Vc;

h2=imagesc(xplot,yplot,Vd);
colormap jet;a = colorbar; daspect([1 1 1])
caxis([-1 1]);axis xy
xlim([100 200]);ylim([-100 150])
set(h2, 'AlphaData' , 1-isnan(Vd))
hold on
title({'Drifters'})
xlabel('x [m]');ylabel('y [m]')
a.Label.String = 'Velocidad [m/s]';
%Curvas de nivel
levels=[1 2 3 4];
[Cc,hc] = contour(bathy.X,bathy.Y,bathy.h,levels);
clabel(Cc,hc)

subplot(1,3,3)

Verr = nanmean(abs(Vd(:)-VG(:)));
Vstd = nanstd(abs(Vd(:)-VG(:)));
DiffV_p = abs(((VG-Vd)./Vd)*100); %Diferencia Porcentual
DiffV = abs(VG-Vd); %Diferencia


h1=imagesc(xplot,yplot,DiffV_p);
colormap jet;a = colorbar; daspect([1 1 1])
caxis([0 200]);axis xy
xlim([100 200]);ylim([-100 150])
set(h1, 'AlphaData' , 1-isnan(LDrift.V_G))
hold on
%quiver(XG(:),YG(:),DIR1,DIR2,0.5,'k')
title(['Mean: ' num2str(round(Verr,3)) ' Std: ' num2str(round(Vstd,3))])
xlabel('x [m]');ylabel('y [m]')
a.Label.String = 'Diferencia Porcentual[%]';

sgtitle(['Velocidad Vy'])

DiffV=abs(Vd(:)-VG(:));
Comp.V.diffVmean=Verr;
Comp.V.diffVstd=Vstd;
Comp.V.diffVrms=sqrt(mean(DiffV .* DiffV, 'omitnan'));

figura2.Position = [0 0 1400 1920]
saveas(figura2,['./Outputs/7-Comparacion/CompV.png'])

%% Correlación entre mediciones in situ y mediciones ópticas
Vi = Vd(:);      % Medición in situ (Drifters) 
Vs = VG(:);     % Medición óptica (Chickadel-Continuidad)

Vi(isnan(Vs)) = nan;
Vs(isnan(Vi)) = nan;
Vi(isnan(Vi)) = [];
Vs(isnan(Vs)) = [];

%Regresión lineal: Yv = beta0 + beta1 * Xv

Xv = [ones(length(Vi),1) Vi];
b = Xv\Vs;
Yv = Xv*b;

x_line=[-1.2:0.1:1.5]';
X_line = [ones(length(x_line),1) x_line];
Y_line = X_line*b;

R2v = 1 - sum((Vs - Yv).^2)/sum((Vs - mean(Vs)).^2)
Comp.V.R2=R2v;

figura5 = figure(5);
clf
scatter(Vi,Vs)
hold on
pv = plot(x_line,Y_line)
pv(1).LineWidth = 1.5;
xlim([-1 1]);ylim([-1 1])

xlabel('Vi [m/s] (mediciones drifters)')
ylabel('Vs [m/s] (método óptico)')

legend({'Data','Regresión Lineal'})
daspect([1 1 1])

txt = [{['Vs = ' num2str(round(b(2),3)) 'Vi + ' num2str(round(b(1),3))]},{['R^2 = ' num2str(round(R2v,3))]}];
text(0.5,0.1,txt)

figura5.Position = [0 0 1400 1920]
saveas(figura5,['./Outputs/7-Comparacion/RegresionV.png'])

%%
figura3 = figure(3);
clf

subplot(1,3,1)
h1=imagesc(xplot,yplot,UG);
colormap jet;a = colorbar; daspect([1 1 1])
caxis([-1 1]);axis xy
xlim([100 200]);ylim([-100 150])
set(h1, 'AlphaData' , 1-isnan(LDrift.V_G))
title({'Chickadel - Continuidad'})
xlabel('x [m]');ylabel('y [m]')
a.Label.String = 'Velocidad [m/s]';
hold on
%Curvas de nivel
levels=[1 2 3 4];
[Cc,hc] = contour(bathy.X,bathy.Y,bathy.h,levels);
clabel(Cc,hc)

subplot(1,3,2)

% Ud=reshape(LDrift.U,flip(size(XG)))';
Ud=LDrift.Uc;

h2=imagesc(xplot,yplot,Ud);
colormap jet;a = colorbar; daspect([1 1 1])
caxis([-1 1]);axis xy
xlim([100 200]);ylim([-100 150])
set(h2, 'AlphaData' , 1-isnan(LDrift.V_G))
hold on
title({'Drifters'})
xlabel('x [m]');ylabel('y [m]')
a.Label.String = 'Velocidad [m/s]';
%Curvas de nivel
levels=[1 2 3 4];
[Cc,hc] = contour(bathy.X,bathy.Y,bathy.h,levels);
clabel(Cc,hc)

subplot(1,3,3)

Uerr=nanmean(abs(Ud(:)-UG(:)));
Ustd = nanstd(abs(Ud(:)-UG(:)));
DiffU_p = abs(((UG-Ud)./Ud)*100); %Diferencia Porcentual
DiffU = abs(UG-Ud); %Diferencia


h1=imagesc(xplot,yplot,DiffU_p);
colormap jet;a = colorbar; daspect([1 1 1])
caxis([0 200]);axis xy
xlim([100 200]);ylim([-100 150])
set(h1, 'AlphaData' , 1-isnan(LDrift.V_G))
hold on
%quiver(XG(:),YG(:),DIR1,DIR2,0.5,'k')
title(['Mean: ' num2str(round(Uerr,3)) ' Std: ' num2str(round(Ustd,3))])
xlabel('x [m]');ylabel('y [m]')
a.Label.String = 'Diferencia Porcentual[%]';

sgtitle(['Velocidad Ux'])


DiffU=abs(Ud(:)-UG(:));
Comp.U.diffUmean=Uerr;
Comp.U.diffUstd=Ustd;
Comp.U.diffUrms=sqrt(mean(DiffU .* DiffU, 'omitnan'));

figura3.Position = [0 0 1400 1920]
saveas(figura3,['./Outputs/7-Comparacion/CompU.png'])

%% Correlación entre mediciones in situ y mediciones ópticas
Ui = Ud(:);      % Medición in situ (Drifters) 
Us = UG(:);     % Medición óptica (Chickadel-Continuidad)

Ui(isnan(Us)) = nan;
Us(isnan(Ui)) = nan;
Ui(isnan(Ui)) = [];
Us(isnan(Us)) = [];

%Regresión lineal: Yv = beta0 + beta1 * Xv

Xu = [ones(length(Ui),1) Ui];
b = Xu\Us;
Yu = Xu*b;

x_line=[-1.2:0.1:1.7]';
X_line = [ones(length(x_line),1) x_line];
Y_line = X_line*b;

R2 = 1 - sum((Us - Yu).^2)/sum((Us - mean(Us)).^2);
Comp.U.R2=R2;

figura6 = figure(6);
clf
scatter(Ui,Us)
hold on
pv = plot(x_line,Y_line)
pv(1).LineWidth = 1.5;
xlim([-1.2 1.5]);ylim([-1.2 1.5])

xlabel('Ui [m/s] (mediciones drifters)')
ylabel('Us [m/s] (método óptico)')

legend({'Data','Regresión Lineal'})
daspect([1 1 1])

txt = [{['Us = ' num2str(round(b(2),3)) 'Ui + ' num2str(round(b(1),3))]},{['R^2 = ' num2str(round(R2,3))]}];
text(0.7,0.4,txt)

figura6.Position = [0 0 1400 1920]
saveas(figura6,['./Outputs/7-Comparacion/RegresionU.png'])

%% Comparacion Direccion

dirD(:,1)=Ui; dirD(:,2)=Vi;
dirD=dirD./sqrt(Ui.^2+Vi.^2);
angD=angle(complex(dirD(:,1),dirD(:,2)))*180/pi();

dirC(:,1)=Us; dirC(:,2)=Vs;
dirC=dirC./sqrt(Us.^2+Vs.^2);
angC=angle(complex(dirC(:,1),dirC(:,2)))*180/pi();

normDeg = mod(angD-angC,360);
DiffAng = min(360-normDeg, normDeg);

diffAmean = mean(DiffAng)
diffAstd  = std(DiffAng)
diffArms  = sqrt(mean(DiffAng .* DiffAng))

Comp.Ang.diffAmean = diffAmean;
Comp.Ang.diffAstd  = diffAstd;
Comp.Ang.diffArms  = diffArms;

%% Save Section

% Output Name

oname=['Coef_Estadisticos'];

% OutPut Directory
odir='./Outputs/7-Comparacion';

% Save file
save([odir '/' oname ],'Comp')



%% GRAFICA DRIFTERS  MAGNITUD-DIRECCION / COMPONENTE U / COMPONENTE V

% figure(7)
% clf
% 
% subplot(1,3,1)
% h2=imagesc(xplot,yplot,LDrift.V_G);
% colormap jet;b = colorbar; daspect([1 1 1])
% caxis([0 1]);axis xy
% xlim([100 200]);ylim([-100 150])
% set(h2, 'AlphaData' , 1-isnan(LDrift.V_G))
% hold on
% quiver(LDrift.Xplot(:),LDrift.Yplot(:),LDrift.DIR1(:),LDrift.DIR2(:),0.3,'k')
% 
% title({'Magnitud Velocidad y Drirección'})
% xlabel('x [m]');ylabel('y [m]')
% b.Label.String = 'Velocidad [m/s]';
% %Curvas de nivel
% levels=[1 2 3 4];
% [Cc,hc] = contour(bathy.X,bathy.Y,bathy.h,levels)
% clabel(Cc,hc)
% 
% 
% subplot(1,3,2)
% h2=imagesc(xplot,yplot,Ud);
% colormap jet;a = colorbar; daspect([1 1 1])
% caxis([-1 1]);axis xy
% xlim([100 200]);ylim([-100 150])
% set(h2, 'AlphaData' , 1-isnan(LDrift.V_G))
% hold on
% title({'Velocidad "u"'})
% xlabel('x [m]');ylabel('y [m]')
% a.Label.String = 'Velocidad [m/s]';
% %Curvas de nivel
% levels=[1 2 3 4];
% [Cc,hc] = contour(bathy.X,bathy.Y,bathy.h,levels)
% clabel(Cc,hc)
% 
% 
% 
% 
% subplot(1,3,3)
% h2=imagesc(xplot,yplot,Vd);
% colormap jet;a = colorbar; daspect([1 1 1])
% caxis([-1 1]);axis xy
% xlim([100 200]);ylim([-100 150])
% set(h2, 'AlphaData' , 1-isnan(Vd))
% hold on
% title({'Velocidad "v"'})
% xlabel('x [m]');ylabel('y [m]')
% a.Label.String = 'Velocidad [m/s]';
% %Curvas de nivel
% levels=[1 2 3 4];
% [Cc,hc] = contour(bathy.X,bathy.Y,bathy.h,levels)
% clabel(Cc,hc)
% 




















%Grafica Direcciones 

% figure(4)
% clf
% 
% subplot(1,3,1)
% h1=imagesc(xplot,yplot,Dg1);
% colormap jet;a = colorbar; daspect([1 1 1])
% caxis([-1 1]);axis xy
% xlim([70 300]);ylim([-300 300])
% set(h1, 'AlphaData' , 1-isnan(Dg1))
% title({'Chickadel - Continuidad'})
% xlabel('x [m]');ylabel('y [m]')
% a.Label.String = 'Velocidad [m/s]';
% %hold on
% 
% subplot(1,3,2)
% 
% D1d=reshape(LDrift.DIR(1,:),flip(size(XG)))';
% 
% h2=imagesc(xplot,yplot,D1d);
% colormap jet;a = colorbar; daspect([1 1 1])
% caxis([-1 1]);axis xy
% xlim([70 300]);ylim([-300 300])
% set(h2, 'AlphaData' , 1-isnan(LDrift.V_G))
% hold on
% title({'Drifters'})
% xlabel('x [m]');ylabel('y [m]')
% a.Label.String = 'Velocidad [m/s]';
% 
% 
% subplot(1,3,3)
% 
% D1err=nanmean(nanmean(abs(D1d-Dg1)));
% 
% h1=imagesc(xplot,yplot,abs(D1d-Dg1));
% colormap jet;a = colorbar; daspect([1 1 1])
% caxis([0 1]);axis xy
% xlim([70 300]);ylim([-300 300])
% set(h1, 'AlphaData' , 1-isnan(LDrift.V_G))
% hold on
% %quiver(XG(:),YG(:),DIR1,DIR2,0.5,'k')
% title(['Diferencia: ' num2str(round(D1err,3))])
% xlabel('x [m]');ylabel('y [m]')
% a.Label.String = 'Velocidad [m/s]';
% 
% sgtitle(['Dirección X'])
% 
% %%
% 
% figure(5)
% clf
% 
% subplot(1,3,1)
% h1=imagesc(xplot,yplot,Dg2);
% colormap jet;a = colorbar; daspect([1 1 1])
% caxis([-1 1]);axis xy
% xlim([70 300]);ylim([-300 300])
% set(h1, 'AlphaData' , 1-isnan(Dg2))
% title({'Chickadel - Continuidad'})
% xlabel('x [m]');ylabel('y [m]')
% a.Label.String = 'Velocidad [m/s]';
% %hold on
% 
% subplot(1,3,2)
% 
% D2d=reshape(LDrift.DIR(2,:),flip(size(XG)))';
% 
% h2=imagesc(xplot,yplot,D2d);
% colormap jet;a = colorbar; daspect([1 1 1])
% caxis([-1 1]);axis xy
% xlim([70 300]);ylim([-300 300])
% set(h2, 'AlphaData' , 1-isnan(LDrift.V_G))
% hold on
% title({'Drifters'})
% xlabel('x [m]');ylabel('y [m]')
% a.Label.String = 'Velocidad [m/s]';
% 
% 
% subplot(1,3,3)
% 
% D2err=nanmean(nanmean(abs(D2d-Dg2)));
% 
% h1=imagesc(xplot,yplot,abs(D2d-Dg2));
% colormap jet;a = colorbar; daspect([1 1 1])
% caxis([0 1]);axis xy
% xlim([70 300]);ylim([-300 300])
% set(h1, 'AlphaData' , 1-isnan(LDrift.V_G))
% hold on
% %quiver(XG(:),YG(:),DIR1,DIR2,0.5,'k')
% title(['Diferencia: ' num2str(round(D2err,3))])
% xlabel('x [m]');ylabel('y [m]')
% a.Label.String = 'Velocidad [m/s]';
% 
% sgtitle(['Dirección Y'])
% 
