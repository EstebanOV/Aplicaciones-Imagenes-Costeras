% ESTIMACIÓN DE LA BATIMETRÍA CON CBATHY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

addpath(genpath('.')) % adds folder dependencies 
% find the path to the cBathy Toolbox
cBathyToolboxPath = what('cBathy-Toolbox');
% cd to the cBathy Toolbox directory to run demo
cd(cBathyToolboxPath.path)

%% Elija el archivo cBathy que desea analizar

load('15V3_Archivo_cBathy-5m.mat')
t=stack.dn.*(24*3600);
data=stack.data;
xyz=stack.xyzAll;
clear stack
save cBathyArrayActual

%% Configuración de INPUTS

edit argus02a

%% RUN the cBathyDemo


clear
tic
% Do a single cBathy analysis, returning the bathy structure (which you
% would normally save somewhere) and displaying the results.  
stationStr = 'argus02a';
stackName = 'cBathyArrayActual';
bathy = analyzeSingleBathyRunNotCIL(stackName, stationStr);
bathy.params.debug.production=1;
toc
%%

plotBathyCollect(bathy)

%% Seccion: Output/Saving

% Output Name
oname=['15V3-bathy-5m-L3'];

% OutPut Directory
odir=['./Outputs/3-Batimetrias/'];

% Save file
save([odir oname],'bathy')

