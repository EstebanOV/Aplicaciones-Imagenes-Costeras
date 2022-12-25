%% GRÁFICA DE RESULTADOS VIDEO_CURRENT_GEN
%==========================================

clc
clear
close all
run Configuracion

%INPUTS
Tw=Tw; % Twindow (segundos)
Ts=Ts; % Tstep   (segundos)

% Archivo de contiene resultados de algoritmo Chickadel
DirOutChick=['./Outputs/2-OutputChickadel/' NameV '_OutputChickadel_Tw' num2str(Tw) '_Ts' num2str(Ts) '.mat'];

load(DirOutChick); 

% Archivo con información de los Timestack
DirTimestack=['./Outputs/1-Instrumentos/' NameV '_Timestack_Ly' num2str(dy) '_dx' num2str(dx) '.mat'];
load(DirTimestack) %Archivo que contiene 3 estructuras:
                                      %InfoTimestack: Coordenadas por pixel y tiempo asociado a cada instrumto
                                      
                                      %TIMESTACK: Esto ya fue utilizado, no debería requerirse nuevamente

% Información del instrumento   
DirIns=['./Outputs/1-Instrumentos/' NameV '_Instrumento_Ly' num2str(dy) '_dx' num2str(dx) '.mat']; 
LoadIns=load(DirIns);
Instr=LoadIns.instr;

%Coordenadas y vector de tiempo asociado a cada frame rectificado 
DirDataF= ['./Inputs/DataFrames_' NameV '.mat'];
LoadData=load(DirDataF);
DataFrames=LoadData.DataFrames;


% Ubicación imágenes rectificadas 
% (utilizado sólo para generar gráficos finales)
imageDirectory= ['C:\Memoria2020\Resultados\Rectificados\' NameV '\MAT'];

% Lista con nombre y dirección de cada imagen
L=string(ls(imageDirectory));
L=L(3:end);
Nframes=length(L);    
for i=1:Nframes
    L(i)=strcat('FrameRec_',num2str(i),'.mat');
end


%Tiempos
TCk=videoCurrentOut(200).InstVbar.t*24*3600; %Vector de tiempo Chikadel tiempo total OBS: Para todos las ubicaciones es el mismo vector
TCk=TCk-TCk(1);                            %tiempo inicializado en primer frame

dTstep=mean(diff(TCk)); %Ventana de tiempo, donde se calcula la velocidad media
TCk2=TCk+dTstep;        %Vector que contiene las ventanas de tiempos acumuladas
dNT=dTstep/0.5;

%% Matriz de aprobación 

%Contrucción de matriz de aprobación: Se acepta o rechaza para cada tiempo
%y ubicación
x=Instr(:,1);
NP=length(x); %Número de puntos
NT=length(videoCurrentOut(200).InstVbar.meanV); %300 %Número de puntos temporales defindos en algoritmo Chikadel 

MAprob=zeros(NP,NT);

for i=1:NP
    for j=1:NT
        sumaCA=sum(CritAcep(i).CA(j,:));
        if sumaCA==3     % Se acepta
            MAprob(i,j)=1;
        else
            MAprob(i,j)=NaN;
        end
    end
end

%% CONSTRUCCIÓN GRILLA X-Y
%-----------------------------

x=Instr(:,1);
y=(Instr(:,2)+Instr(:,3))*0.5;
dx=LoadIns.dx; % [m]  Input
dy=LoadIns.Ly; % [m]  Input

Npx=(max(x)-min(x))/dx; %Numero de Particiones en x
Npy=(max(y)-min(y))/dy; %Numero de Particiones en y
xlin=linspace(min(x),max(x),Npx+1);
ylin=linspace(min(y),max(y),Npy+1);
[X,Y]=meshgrid(xlin,ylin); % Grilla X e Y
Yf=fliplr(Y')';

NP=length(x); %Número de puntos
NT=length(videoCurrentOut(200).InstVbar.meanV); %Número de puntos temporales defindos en algoritmo Chikadel
Mv=zeros(NP,NT);

for i=1:NP
    for j=1:NT
    if MAprob(i,j)==1
        Mv(i,j)=videoCurrentOut(i).InstVbar.meanV(j);       
    else
        Mv(i,j)=NaN;
    end
    end
end


% CONSTRUCCIÓN GRILLA VELOCIDAD 
%--------------------------------

for i=1:NT
VelGrilla(i).MV=(1)*griddata(x,y,Mv(:,i),X,Y,'linear'); % AACÁAAAAAAA PUUUUUUUSSSEEEEEEEEEEEEEEEEE EEEEEEEEEEEEELLLLLLLLLLLLL -------11111111111111
                                                         %ARREGLAR CUANDO
                                                         %EJECUTE EL CODIGO
                                                         %DESDE EL
                                                         %PRINCIPIO                                                        
end


%% CONTRUCCIÓN FRAMES PROMEDIADOS SEGÚN VENTANA DE TIEMPO (Tstep)
%----------------------------------------------------------------
%Frames promediados quedan almacenados en Variable Iprom

% c1=0;
% for i=1:length(TCk)
%     iTimex=double(load(strcat(imageDirectory, L(100))).Irect).*0;
%     c1=1+dNT*(i-1);
%     c2=dNT*i;
%     for j=c1:c2
%         iTimex=iTimex+double(load(strcat(imageDirectory, L(j))).Irect);
%         if j==c2
%         IProm{i}=uint8(iTimex./dNT);   %Cálculo promedio 
%         end
%     end
% end

%%
% GRÁFICA RESULTADOS: VELOCIDAD MEDIA Y FRAMES PROMEDIADOS 
%----------------------------------------------------------

% figura1=figure(1)
% U=zeros(length(ylin),length(xlin));
% 
% for i=1:NT
% 
% %GRÁFICA MATRIZ DE VELOCIDAD 
% %---------------------------
% 
% subplot(1,2,2)
% 
% imagesc(xlin,ylin,VelGrilla(i).MV)
% a = colorbar;
% caxis([-3 3])
% axis xy
% xlim([100 360])
% ylim([-270 260])
% title(['Vy , T = '+string(TCk2(i))+'s'])
% xlabel('Cross-shore position [m]')
% ylabel('Alongshore position[m]')
% a.Label.String = 'Velocidad [m/s]';
% daspect([1 1 1])
% 
% hold on;
% 
% %Grafica Velocidad como Campo Vectorial
% quiver(X,Y,U,VelGrilla(i).MV,'r','LineWidth',1,'AutoScaleFactor',1.5)
% xlim([115 360])
% ylim([-270 250])
% 
% 
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
% xlim([100 370])
% ylim([-270 260])
% title(['Imágenes Promediadas: \Deltat = ' num2str(dTstep) 's'])
% 
% 
% 
% sgtitle(['Vy , Tw = ' num2str(Tw) ' (' num2str(Tw*0.5) 's) , Ts = ' num2str(Ts) ' (' num2str(Ts*0.5) 's)'])
% %saveas(figura1,['./Resultados/PlotVy_Tw' num2str(Tw) '_Ts' num2str(Ts) '/Vy_T' num2str(TCk2(i)) 's.png'])
% 
% pause(0.1)
% end

%% VELOCIDAD PROMEDIADA
% -----------------------------

%Promedio y desviación estandar
dim=size(VelGrilla(1).MV);
Ni=dim(2);
Nj=dim(1);

for i=1:Ni
    for j=1:Nj
        for k =1:NT
            ListaVij(k)=VelGrilla(k).MV(j,i);
            
        end
        
        GCorrientes.Vy_mean(j,i)=mean(ListaVij,'omitnan');
        GCorrientes.Vy_st(j,i)=std(ListaVij,'omitnan');
        
    end
end

%% GRÁFICA VELOCIDAD PROMEDIO VY

figura2=figure(2)
U=zeros(length(ylin),length(xlin));

%GRÁFICA MATRIZ DE VELOCIDAD 
%---------------------------

%Matriz Velocidad
subplot(1,2,1)

imagesc(xlin,ylin,GCorrientes.Vy_mean)
a = colorbar;
caxis([-1.5 1.5])
axis xy
xlim([65 360])
ylim([-270 260])
title(['Mean Vy , Tw = ' num2str(Tw)])
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Velocidad [m/s]';
daspect([1 1 1])

hold on;


%Matriz Desviación Estandar
subplot(1,2,2)

imagesc(xlin,ylin,GCorrientes.Vy_st)
a = colorbar;
caxis([0 1])
axis xy
xlim([65 360])
ylim([-270 260])
title(['Std Vy , Tw = ' num2str(Tw)])
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Velocidad [m/s]';
daspect([1 1 1])

sgtitle(['Mean Vy , Tw = ' num2str(Tw) ' s , Ts = ' num2str(Ts) ' s'])
figura2.Position = [0 0 1400 1920]
saveas(figura2,['./Outputs/3-Velocidades/' NameV 'PlotVy_Tw' num2str(Tw) '_Ts' num2str(Ts) '.png'])


%% BINARIZACIÓN Y FILTRADO DE IMAGEN


dim=size(GCorrientes.Vy_mean);
Vy_mean=GCorrientes.Vy_mean;
Ni=dim(2);
Nj=dim(1);
for i=1:Ni
   for j=1:Nj
       if isnan(Vy_mean(j,i))   
           MV_filt(j,i)=0;
       else
          MV_filt(j,i)=1;
       end
   end
end


I1=mat2gray(MV_filt); %Imagen a escala de grises
%I1=GCorrientes.Vy_mean
I2=imbinarize(I1);                %Imagen binarizada
%I2=imbinarize(I1,'adaptive','ForegroundPolarity','dark','Sensitivity',0.4);
I3=imcomplement(I2);              %Imagen Complementaria
I4=imfill(I3,4,'holes');          %Relleno de Regiones
I4=imcomplement(I4); 
I4 = bwareaopen(I4, 4);           %Quitar objetos pequeños de la imagen
I5=imcomplement(I4);              %Imagen Complementaria
I5 = bwareaopen(I5, 5);           %Quitar objetos pequeños de la imagen
I5 = imcomplement(I5);            %Imagen Complementaria
se6 = strel('disk',1);
% se7 = strel('disk',1);
% I6=imclose(I4,se6);
I6=imcomplement(I5);
I8=imopen(I6,se6);               %Apertura morfologica de la imagen
I8=imcomplement(I8);

%Con filtrado Gausiano
Num=20;
IR=imresize(I8, Num); %I5 en vez de I8 (otra opción)
I9=imguidedfilter(IR,'NeighborhoodSize',[2 2]);
I9=imresize(I9, 1/Num);
I10=imguidedfilter(IR,'DegreeOfSmoothing',700);
I10=imresize(I10, 1/Num);

%Filtrado Gaussiano Complementario
I10c = imcomplement(I10);
IR=imresize(I10c, Num);
I99=imguidedfilter(IR,'NeighborhoodSize',[2 2]);
I99=imresize(I99, 1/Num);
I100=imguidedfilter(IR,'DegreeOfSmoothing',700);
I100=imresize(I100, 1/Num);
I100 = imcomplement(I100);
I100=imopen(I100,se6);

I100 = imcomplement(I100);
I100 = bwareaopen(I100, 10); 
I100 = imcomplement(I100);

% Opcional (Comente o descomente en caso de estar conforme con los
% resultados)
I100=imcomplement(I100);
I100=imopen(I100,se6);               %Apertura morfologica de la imagen
I100=imcomplement(I100);

%Mas filtros
% Num2=3;
% I11=imresize(I10,Num2);
% I11=medfilt2(I11,[2 2])
% I11=imresize(I11,1/Num2);
% I11=imresize(I11,Num2);
% I11=medfilt2(I11,[2 2])
% I11=imresize(I11,1/Num2);
% I11=imcomplement(I11);
% I11 = bwareaopen(I11, 2);
% I11=imresize(I11,Num2);
% I11=medfilt2(I11,[2 2])
% I11=imresize(I11,1/Num2);

figure3=figure(3)

subplot(1,3,1)

% imagesc(xlin,ylin,GCorrientes.Vy_mean)
% a = colorbar;
% caxis([-1.5 1.5])
% axis xy
% xlim([65 360])
% ylim([-270 260])
% title(['Mean Vy , Tw = ' num2str(Tw)])
% xlabel('Cross-shore position [m]')
% ylabel('Alongshore position[m]')
% a.Label.String = 'Velocidad [m/s]';
% daspect([1 1 1])

imshow(flipud(I1))

subplot(1,3,2)
%imshow(flipud(I9))
imshow(flipud(I5))

subplot(1,3,3)
imshow(flipud(I100))

sgtitle(['Binarización, Operaciones Morfológicas y Filtrado'])
figure3.Position = [0 0 1400 1920]
saveas(figure3,['./Outputs/3-Velocidades/' NameV '_PlotVy_Filtro_Tw' num2str(Tw) '_Ts' num2str(Ts) '.png'])

%% %% Interpolación

% Lista X,Y,Velocidad
LV=reshape(GCorrientes.Vy_mean,[],1);
LVV=reshape(GCorrientes.Vy_mean,[],1);
Lx=LoadIns.instr(:,1);
Ly=LoadIns.instr(:,2);

% Eliminar valores Nan
for i=1:length(LV)
   if isnan(LV(i))
       Lx(i)=NaN;
       Ly(i)=NaN;
   end
end

LV=LV(~isnan(LV));
Lx=Lx(~isnan(Lx));
Ly=Ly(~isnan(Ly));


%Interpolación
%[x_grilla,y_grilla]=meshgrid(x_g1,y_g1);
x=LoadIns.instr(:,1);
y=LoadIns.instr(:,2);
%MV_int=griddata(Lx,Ly,LV,x,y,'linear'); %Interpolación Lineal
FV=scatteredInterpolant(Lx,Ly,LV);
FV.Method='linear';
MV_int=FV(x,y);

sizze=size(GCorrientes.Vy_mean);
Nfil=sizze(1);
Ncol=sizze(2);
MV_int=reshape(MV_int,[Nfil,Ncol]);

% Filtro 

for i=1:Ni
   for j=1:Nj
       if I100(j,i)   
           MV_filt(j,i)=MV_int(j,i);
       else
          MV_filt(j,i)=NaN;
       end
   end
end

GCorrientes.VyF_mean=MV_filt;
% Grafica

figure4=figure(4)

subplot(1,2,1)

imagesc(xlin,ylin,GCorrientes.Vy_mean)
a = colorbar;
caxis([-1.5 1.5])
axis xy
xlim([65 360])
ylim([-270 260])
title(['Mean Vy , Tw = ' num2str(Tw)])
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Velocidad [m/s]';
daspect([1 1 1])


subplot(1,2,2)
imagesc(xlin,ylin,MV_filt)
a = colorbar;
caxis([-1.5 1.5])
axis xy
xlim([65 360])
ylim([-270 260])
title(['Mean Vy , Tw = ' num2str(Tw) ' (FILTRADO)'])
xlabel('Cross-shore position [m]')
ylabel('Alongshore position[m]')
a.Label.String = 'Velocidad [m/s]';
daspect([1 1 1])

figure4.Position = [0 0 1400 1920]
saveas(figure4,['./Outputs/3-Velocidades/' NameV '_PlotVy_scF_Tw' num2str(Tw) '_Ts' num2str(Ts) '.png'])





%% Sección: Output/Saving

COutVy.GrillaX=X;
COutVy.GrillaY=Y;
COutVy.x_list=xlin;
COutVy.y_list=ylin;
COutVy.dx=dx;
COutVy.dy=dy;
COutVy.TCk=TCk;
COutVy.GrillaVy=VelGrilla;
COutVy.GrillaVyF=MV_filt;
%COutVy.FramesMean=IProm; %Descomentar si quiere guardar los frames
%promediados

% Output Name
oname=[NameV 'C_Output_Vydx' num2str(dx) 'dy' num2str(dy) '_Tw' num2str(Tw) 'Ts' num2str(Ts)];

% OutPut Directory
odir='./Outputs/3-Velocidades';

% Save file
save([odir '/' oname ],'COutVy','GCorrientes')

close all


