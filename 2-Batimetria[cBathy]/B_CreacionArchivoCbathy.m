%% B) Creacion Archivo cBathy


clear
addpath(genpath('.'))

%% Seccion 1: Inputs (Pixel Instruments)

load('./Outputs/1-PixelInstruments/15V3_5m_pixInst.mat')
load('./Inputs/DataFrames_15V3.mat')

%% Seccion 2: Main Code

%Tipo de Instrumento

stack.inst.type='matrix';
stack.inst.xyz=[];
stack.inst.name='cBathyArray';
stack.inst.ShortName='mBW';
stack.inst.x=[pixInst.xlim(1) pixInst.dx pixInst.xlim(2)];
stack.inst.y=[pixInst.ylim(1) pixInst.dy pixInst.ylim(2)];
stack.inst.z=0;

NxNy=size(pixInst.X);
Nx=NxNy(2);
Ny=NxNy(1);
count=0;
for i=1:Nx
   for j=1:Ny
       count=count+1;
       xyzAll(count,:)=[pixInst.X(j,i) pixInst.Y(j,i) 0];
       
   end
end

stack.inst.xyzAll=xyzAll;
deltaT=0.5/(3600*24); % segundos

NT=size(pixInst.Igray);
NT=NT(3);
TiempoF=(NT-1)*deltaT;
dn=[0:deltaT:TiempoF];
dn=dn+DataFrames.Time(1);
stack.dn=dn;
stack.xyzAll=xyzAll;

for k=1:NT
   II=pixInst.Igray(:,:,k);
   count=0;
   for i=1:Nx
       for j=1:Ny
           count=count+1;
           data(k,count)=II(j,i);
       end      
   end
end

stack.data=data;

%% Seccion 3: Output/Saving

% Output Name
oname=['15V3_Archivo_cBathy-5m'];

% OutPut Directory
odir=['./Outputs/2-Archivos-cBathy/'];

% Save file
save([odir oname],'stack')