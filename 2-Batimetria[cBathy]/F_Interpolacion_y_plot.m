%% GRAFICA BATIMETRÍA E INTERPOLACIÓN

clear
%Profundidad filtro Kalman
Bf = load('./BathysSmoothed/BathysSmoothed_6/15V3-bathy-5m-L3.mat');

%Profundidad 
h = Bf.bathy.runningAverage.h;

%Error
hErr = Bf.bathy.runningAverage.hErr;

%Remoción según umbral de error
threshErr = 2;
h(hErr>threshErr) = nan;
H = h(:); 

%Coordendas y Grilla Bathy
[X,Y] = meshgrid(Bf.bathy.xm,Bf.bathy.ym);
Xi1 = X(:); Yi = Y(:);

%Coordenadas y Grillas para interpolar ( limites x = [80 - 300] )
Xi=Xi1(Xi1<301); Yi=Yi(Xi1<301); H=H(Xi1<301); 
xg = Bf.bathy.xm; xg = xg(xg<301); yg = Bf.bathy.ym;
[Xg,Yg] = meshgrid(xg,yg);

%Interpolación
H=reshape(H,size(Xg));
hi = inpaint_nans(H); %Propfundidad interpolada

%% Graficas

figure(1)

%Grafica profundidad cBathy (con filtro de Kalman)
subplot(1,2,1)
imagesc(Bf.bathy.xm,Bf.bathy.ym,-h); grid on
colormap parula
a = colorbar; a.Label.String = 'Profundidad [m]';
caxis([-10 0])
axis xy
xlim([80 300]); ylim([-300 300])
axis xy; xlabel('x (m)'); ylabel('y (m)');
title(['Profundidad cBathy'])
daspect([1 1 1])
hold on
contour(X,Y,h,'ShowText','on')


%Grafica profundidad cBathy con interpolación
subplot(1,2,2)
imagesc(xg,yg,-hi); grid on
colormap parula; a = colorbar;
caxis([-10 0])
axis xy
xlim([80 300]) ;ylim([-300 300])
axis xy; xlabel('x (m)'); ylabel('y (m)');
title(['Prof cBathy con Interpolación'])
a.Label.String = 'Profundidad [m]';
daspect([1 1 1])
hold on
contour(Xg,Yg,hi,'ShowText','on')

%% Seccion: Output/Saving

%Variables a Guardar

bathy.h=hi;
bathy.X=Xg;
bathy.Y=Yg;
bathy.xm=xg;
bathy.ym=yg;
% Output Name
oname=['Bathy-5m-Las-Cruces'];

% OutPut Directory
odir=['./Outputs/4-Bathy-Interpolada/'];

% Save file
save([odir oname],'bathy')

























% Papelera de Codigos
%  INTERPOLACIÓN

% H = h(:);  H = H(~isnan(H));
% Xi = X(:); Xi = Xi(~isnan(H));
% Yi = Y(:); Yi = Yi(~isnan(H));
% 
% x = Bf.bathy.xm; x = x(x<301);
% y = Bf.bathy.ym;
% [Xg,Yg] = meshgrid(x,y);
% 
% %Limites x
% Xi=Xi(Xi<301); Yi=Yi(Xi<301); H=H(Xi<301);
% 
% %%
% %Interpolación
% %hf = griddata(Xi,Yi,H,Xg,Yg,'linear'); 
% hf = interp2(Xi,Yi,H,Xg,Yg,'spline'); 
% 
% %Grafica
% subplot(1,2,2)
% imagesc(Bf.bathy.xm,Bf.bathy.ym,-hf); grid on
% colormap parula
% a = colorbar;
% caxis([-10 0])
% axis xy
% xlim([80 300])
% ylim([-300 300])
% axis xy; xlabel('x (m)'); ylabel('y (m)');
% title(['cBathy Filtro Kalman (Int)'])
% a.Label.String = 'Profundidad [m]';
% daspect([1 1 1])
% hold on
% contour(Xg,Yg,hf,'ShowText','on')

%
% 
% F = scatteredInterpolant(Xi,Yi,H); %Función de interpolación
% F.Method = 'linear';
% HF= F(X(:),Y(:));
% 
% hf=reshape(HF,size(h));





