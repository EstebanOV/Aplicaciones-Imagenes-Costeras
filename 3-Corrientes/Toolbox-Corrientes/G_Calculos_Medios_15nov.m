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
B15V1 = load(['./Outputs/2-OutputChickadel/15V1_OutputChickadel_Tw' num2str(Tw) '_Ts' num2str(Ts) '.mat']);
B15V2 = load(['./Outputs/2-OutputChickadel/15V2_OutputChickadel_Tw' num2str(Tw) '_Ts' num2str(Ts) '.mat']);
B15V3 = load(['./Outputs/2-OutputChickadel/15V3_OutputChickadel_Tw' num2str(Tw) '_Ts' num2str(Ts) '.mat']);

%Archivos arreglos velocidad Vy (Outputs Scripts C)
C15V1 = load(['./Outputs/3-Velocidades/15V1C_Output_Vydx' num2str(dx) 'dy' num2str(dy) '_Tw' num2str(Tw) 'Ts' num2str(Ts) '.mat']);
C15V2 = load(['./Outputs/3-Velocidades/15V2C_Output_Vydx' num2str(dx) 'dy' num2str(dy) '_Tw' num2str(Tw) 'Ts' num2str(Ts) '.mat']);
C15V3 = load(['./Outputs/3-Velocidades/15V3C_Output_Vydx' num2str(dx) 'dy' num2str(dy) '_Tw' num2str(Tw) 'Ts' num2str(Ts) '.mat']);

%Cargar Batimetría
load('./Inputs/Bathy-5m-Las-Cruces.mat')

%Velocidades
V5=C15V1.GCorrientes.Vy_mean;
V6=C15V2.GCorrientes.Vy_mean;
V7=C15V3.GCorrientes.Vy_mean;


%Coordenadas
x_g1 = C15V1.COutVy.x_list; 
y_g1 = C15V1.COutVy.y_list; 

%Grilla Vy
[x_grilla1,y_grilla1]=meshgrid(x_g1,y_g1);

% Información del instrumento
Name='15V1'; %Todos son iguales
DirIns=['./Outputs/1-Instrumentos/' Name '_Instrumento_Ly' num2str(dy) '_dx' num2str(dx) '.mat']; 
LoadIns=load(DirIns);Instr=LoadIns.instr;



% %Filtro Velocidades
% V1(V1>1.5)=1.5; V1(V1<-1.5)=-1.5;  %Para otros datos puede ser útil

% Promedio Velocidades
VM=reshape(nanmean([V5(:)';V6(:)';V7(:)']),size(V5));
VMstd=reshape(nanstd([V5(:)';V6(:)';V7(:)']),size(V5));

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
title('División estándar')
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Velocidad [m/s]';
set(bbp, 'AlphaData' , 1-isnan(VMstd))

sgtitle(['Velocidad Media (15 nov)'])
figura1.Position = [0 0 1400 1920]
saveas(figura1,['./Outputs/5-Graficas-Promedio7v/VyMean-Std-15nov.png'])
%%
%Filtro de Velocidades Altas
VMn = VM; VMn(abs(VMn)>0.6) = nan; %VMn(VMn>0.4) = nan; 0.56

%Interpolación en datos nan
VMi= inpaint_nans(VMn);

%Filtro Gaussiano
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

sgtitle(['Interpolación y Filtrado (15 nov)'])
figura2.Position = [0 0 1400 1920]
saveas(figura2,['./Outputs/5-Graficas-Promedio7v/Vy-Int_Filt-15nov.png'])

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

% Cálculo Ux

N2=length(x_g2);
M2=length(y_g2);
h=Depth_g2*(-1);       % Profundidad [m]

VMi=VMg;               % Velocidad Vy [m/s]
U=nan(M2,N2);          % Inicialización Matriz U 
U(:,1)=0;              % Condicion de Bord1


for i=1:(N2-1)
   for j=1:(M2-1)
       
       hy1=mean([h(j,i) h(j,i+1)]);
       hy2=mean([h(j+1,i) h(j+1,i+1)]);
       hx1=mean([h(j,i) h(j+1,i)]);
       hx2=mean([h(j,i+1) h(j+1,i+1)]);
       V1=VMi(j,i);
       V2=VMi(j+1,i);
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


Ux=U; Ux(M2,:)=[];
Uxt2=U;
Uxt2(:,1)=[];

%Grilla Ux
xu = x_g2; yu = y_g2;
yu = yu + dy/2 ; yu(length(yu)) = [];
[xu_g,yu_g]=meshgrid(xu,yu);

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
sgtitle(['Campo de Velocidades (15 nov)'])
figura4.Position = [0 0 1400 1920]
saveas(figura4,['./Outputs/5-Graficas-Promedio7v/Vy-Ux-15nov.png'])

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

sgtitle(['Campo de Velocidades (15 nov)'])
figura5.Position = [0 0 1400 1920]
saveas(figura5,['./Outputs/5-Graficas-Promedio7v/Corrientes-15nov.png'])

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
oname=[Name '_CalculosMedios-15nov'];

% OutPut Directory
odir='./Outputs/6-Valores-Medios';

% Save file
save([odir '/' oname ],'VM','UM','xu','yu','Xu','Yu','xv','yv','Xv','Yv','VarPlot')
