%% CREACIÓN DE INSTRUMENTOS

clc
clear
close all
run Configuracion

%% Seccion 1: INPUTS

% Direccion de Carpetas y Archivos

% Direccion de Carpetas y Archivos
DirFrames =['C:\Memoria2020\Resultados\Rectificados\' NameV '\MAT']; %Carpeta donde se ubican los frames
DirDataF= ['./Inputs/DataFrames_' NameV]; %Archivo que contiene información de rectificación
Prefijo='FrameRec_'; %Prefijo del nombre de los frames


% Dimensiones Instrumentos vBar

Ly=dy;          % [m]
dx=dx;          % [m]
Xlim=Xlim;
Ylim=Ylim;
z = 0;

%% Sección 2: CARGA DE ARCHIVOS Y FRAMES

% Nombres Frames

clipFns = dir(DirFrames);
Nf=length(clipFns)-2; %Número de Frames
for i=1:Nf
   ListNf{i}=strcat(Prefijo,num2str(i),'.mat'); %Lista con todos los nombres de los frames
end

% Datos Rectificación

LoadData=load(DirDataF);
DataFrames=LoadData.DataFrames;

 % Coordenda asociada a las imágenes
 xFr=DataFrames.Coord.x;
 yFr=DataFrames.Coord.y;
 yFrFl=fliplr(yFr);


% Carga frame de ejemplo

LoadF1=load(strcat(DirFrames,'\',ListNf{1}));
Frame1=LoadF1.Irect;

%% Seccion 3: DEFINICIÓN INSTRUMENTOS VBAR

x=[Xlim(1):dx:Xlim(2)];
y=[Ylim(1):Ly:Ylim(2)];

Nx=(Xlim(2)-Xlim(1))/dx;
Ny=(Ylim(2)-Ylim(1))/Ly;
cnt = 1;

for i = 1: Nx
    for j = 1:(Ny)
    instr(cnt,1)=x(i);
    instr(cnt,2)=y(j);
    instr(cnt,3)=y(j+1);
    cnt = cnt+1;
    end
end

NIns=length(instr); % Número de instrumentos

%Gráfica de instrumentos

figura10=figure(10); clf
axis xy;axis image
hold on
imagesc(DataFrames.Coord.x(1,:),DataFrames.Coord.y(1,:),uint8(Frame1));
for i = 1: length(instr)
    line([instr(i,1),instr(i,1)], [instr(i,2),instr(i,3)], 'Color', 'r','LineWidth', 2);    
end
plot(instr(:,1),instr(:,2),'.k','MarkerSize',7);
plot(instr(:,1),instr(:,3),'.k','MarkerSize',7);
xlabel('cross-shore position [m]')
ylabel('alongshore position [m]')
title(['Instrumento Malla vBar  Ly=' num2str(Ly) 'm  dx=' num2str(dx) 'm ']);
hold off
figura10.Position = [0 0 1400 1920]

%% Seccion 4: Creación de Instrumentos

% Grafica de ejemplo

 figure(2)
 Frame=double(flipud(rgb2gray(uint8(Frame1))));
 imagesc(Frame)
 colormap(gray)

% Cálculo rellenando Timestack

%h = waitbar(0,'Creando TIMESTACK');

%tic
for i=1:Nf
     LoadFi=load(strcat(DirFrames,'\',ListNf{i}));
     Framei=double(flipud(rgb2gray(uint8(LoadFi.Irect))));
     
     for j=1:NIns
         
         x1=instr(j,1);
         y1=instr(j,2);
         x2=instr(j,1);
         y2=instr(j,3);

         indx1=find(xFr(1,:)==x1);
         indy1=find(yFrFl(1,:)==y1);
         indx2=find(xFr(1,:)==x2);
         indy2=find(yFrFl(1,:)==y2);        
         
         TIMESTACK(j).stacki(i,:)=fliplr(Framei(indy2:indy1,indx1)'); 
          
    end
    
    %waitbar(i/(Nf),h)
end
%waitbar(1,h,'El proceso ha finalizado correctamente')
%toc
%% Sección 5: Poducto Final vBar Instruments


for j=1:NIns
    
    x1=instr(j,1);
    y1=instr(j,2);
    x2=instr(j,1);
    y2=instr(j,3);

    indx1=find(xFr(1,:)==x1);
    indy1=find(yFr(1,:)==y1);
    indx2=find(xFr(1,:)==x2);
    indy2=find(yFr(1,:)==y2);
    
    InfoTimestack(j).xyzAll(:,2)=yFr(1,indy1:indy2)';
    InfoTimestack(j).xyzAll(:,1)=repmat(x1,1,length(InfoTimestack(j).xyzAll(:,2)))';
    InfoTimestack(j).xyzAll(:,3)=repmat(0,1,length(InfoTimestack(j).xyzAll(:,2)))';
    InfoTimestack(j).time=DataFrames.Time;
    InfoTimestack(j).dyFrame=mean(diff(InfoTimestack(j).xyzAll(:,2)));
    
end

%% Sección 6: Output/Saving

% Output Name
oname1=[NameV '_Instrumento_Ly' num2str(Ly) '_dx' num2str(dx)];
oname2=[NameV '_Timestack_Ly' num2str(Ly) '_dx' num2str(dx)]; 

% OutPut Directory
odir='./Outputs/1-Instrumentos';

% Save file
save([odir '/' oname1 ],'instr','Ly','dx','Xlim','Ylim','z','ListNf')
save([odir '/' oname2 ],'TIMESTACK','InfoTimestack')
saveas(figura10,['./Outputs/1-Instrumentos/' NameV '_Ins_vBar_Ly' num2str(Ly) 'm_dx' num2str(dx) 'm.png'])
