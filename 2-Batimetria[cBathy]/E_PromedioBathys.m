%% Promedio Batimetrías

%Load Batimetrias
Bathy_1 = load('./Outputs/3-Batimetrias/12V1-bathy-5m-L3.MAT');
Bathy_2 = load('./Outputs/3-Batimetrias/13V1-bathy-5m-L3.MAT');
Bathy_3 = load('./Outputs/3-Batimetrias/13V2-bathy-5m-L3.MAT');
Bathy_4 = load('./Outputs/3-Batimetrias/13V3-bathy-5m-L3.MAT');
Bathy_5 = load('./Outputs/3-Batimetrias/15V1-bathy-5m-L3.MAT');
Bathy_6 = load('./Outputs/3-Batimetrias/15V2-bathy-5m-L3.MAT');
Bathy_7 = load('./Outputs/3-Batimetrias/15V3-bathy-5m-L3.MAT');

%Profundidades
h1 = Bathy_1.bathy.fCombined.h;
h2 = Bathy_2.bathy.fCombined.h;
h3 = Bathy_3.bathy.fCombined.h;
h4 = Bathy_4.bathy.fCombined.h;
h5 = Bathy_5.bathy.fCombined.h;
h6 = Bathy_6.bathy.fCombined.h;
h7 = Bathy_7.bathy.fCombined.h;

%Error
hErr1 = Bathy_1.bathy.fCombined.hErr;
hErr2 = Bathy_2.bathy.fCombined.hErr;
hErr3 = Bathy_3.bathy.fCombined.hErr;
hErr4 = Bathy_4.bathy.fCombined.hErr;
hErr5 = Bathy_5.bathy.fCombined.hErr;
hErr6 = Bathy_6.bathy.fCombined.hErr;
hErr7 = Bathy_7.bathy.fCombined.hErr;

%% Remoción resultados con un error mayor al umbral especificado

threshErr = 1;

h1(hErr1>threshErr) = nan;
h2(hErr2>threshErr) = nan;
h3(hErr3>threshErr) = nan;
h4(hErr4>threshErr) = nan;
h5(hErr5>threshErr) = nan;
h6(hErr6>threshErr) = nan;
h7(hErr7>threshErr) = nan;

% Promedio 

H1=h1(:)'; H2=h2(:)'; H3=h3(:)'; H4=h4(:)';
H5=h5(:)'; H6=h6(:)'; H7=h7(:)';

Hm = nanmean([H1;H2;H3;H4;H5;H6;H7])'; 
hm = reshape(Hm,size(h1));

%% Grafica profundidades promedios

figure(1)
h=imagesc(Bathy_1.bathy.xm,Bathy_1.bathy.ym,-hm); grid on
colormap parula
a = colorbar;
caxis([-10 0])
axis xy
xlim([80 300])
ylim([-300 300])
axis xy; xlabel('x (m)'); ylabel('y (m)');
title(['cBathy promedio'])
a.Label.String = 'Profundidad [m]';
daspect([1 1 1])
