% Conversión de Archivo .mat a .txt
clear
load('refPOINT_UTM_12.mat')

%% Guardar todo
GCP_UTM_V12=struct2table(gcpGlobalC);
writetable(GCP_UTM_V12)
type 'GCP_UTM_V12.txt'

%% Guardar solo los valores 

% GCP_UTM_V12(:,1)=T.num;
% GCP_UTM_V12(:,2)=T.x;
% GCP_UTM_V12(:,3)=T.y;
% GCP_UTM_V12(:,4)=T.z;
% writematrix(GCP_UTM_V12)

