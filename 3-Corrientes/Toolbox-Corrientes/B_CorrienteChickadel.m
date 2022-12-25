% Description:
% This script applies the videoCurrentGen code to compute optical currents
% for an example video-currents stack structure.
%
% Here, the expected fields of the stack structure are:
% inst: 
%        type: 'line'
%        xyz: [x0 y0 z; x1 y1 z], where (x0,y0,z) and (x1,y1,z) are the
%              endpoints of pixel line, e.g. a 20-m alongshore line:
%              xyz: [115 640 0; 115 660 0], size [3x2]
%        name: 'vBar115'
% dn: time vector (matlab datenum) with size [Mx1s]
% xyzAll: vector with x,y,z coordinates of points along pixel line
%        instrument, size [Nx3]
% data: grayscale image of timestack for one linear pixel instrument,
%        with size [MxN]
% Note on dimensions:
%        M is the length of the whole record, e.g., for a 15 min video,
%           subsampled at 2 Hz, M~1801
%           The timestack will be analyzed in windows Twin (e.g., 64 s)
%           for the video-current estimation in videoCurrentGen
%        N is the length of pixel array given by points (x,y,z), e.g.,
%           for a 20 m pixel array with .2 m spacing, N~101
% Note on instruments:
%        If there are multiple pixel instruments, adapt this script to loop
%        through stackstruct(1), stackstruct(2), etc.
%
% (This should be the same form as the Aerielle video demo output from the
%  2017 bootcamp.)
%% INPUTS
clc
clear
close all
run Configuracion

DirTimestack=['./Outputs/1-Instrumentos/' NameV '_Timestack_Ly' num2str(dy) '_dx' num2str(dx) '.mat']
% Load stack structure
stackStruct=load(DirTimestack);
si=size(stackStruct.TIMESTACK);

%%

ce=0;
for j=1:si(1,2);
% Grab stack: grayscale image of timestack, size [MxN]
stack = stackStruct.TIMESTACK(j).stacki;

% Grab time: time line (starting from zero) of stack, size [1xM]
time = stackStruct.InfoTimestack(j).time;

% Grab xy: x,y position [Nx2] of each pixel in a dimension (equally spaced)
xyz = stackStruct.InfoTimestack(j).xyzAll;
xy = xyz(:,1:2);

% Set Twin: the time length of the FFT window (in points)
% For 2 Hz data, Twin = 128 will yield a 64s average current

Twin = Tw*2; % x2 por los 2 Hz

% Set Tstep: time length to step the window (in points)
% For 2 Hz data, Tstep = 64 will yield a current estimate every 32 s

Tstep = Ts*2; 

% Set vBounds: [minV maxV], units m/s or vector of desired velocity steps,
% Set this to empty, [], to use defaults [-3 3]
vB = [];

% Set fkBounds = [fmin fmax kmin kmax], vector of frequency and wavenumber
%      bounds energy out of side of these bounds will be set to 0.  Useful
%      to eliminate some of the wave contamination that leaks in.  Set this
%      to empty, [], to use defaults.

fkB = [inf 0 km 2]; % sólo se cambia el valor de kmin, los demas son los valores por defecto

% Set {plotFlag}: optional, if true (~=0) will display a running plot of
% %     the data processing
plotFlag = 0;

% RUN Video Current Gen
%j

% if sum(sum(stack))< 2200000
%     videoCurrentOut(j).InstVbar=nan;
% %  elseif   j==9 | j==114 |  j==479 |  j==531 | j==532 | j==584  %13v1?                  %j==107 %Vuelo 15V3  % 2481320   %Vuelo 13V1
% %      videoCurrentOut(j).InstVbar=nan;  
%  else
try
    videoCurrentOut(j).InstVbar = videoCurrentGen(stack, time, xy, vB, fkB,Twin, Tstep, plotFlag);
catch
    videoCurrentOut(j).InstVbar=nan;
    ce=ce+1;      %Contador de errores
    errores(ce)=j;
end
end 

%  OUTPUT fields in videoCurrentOut returned:
%    meanV - video current estimate of mean current for each time window
%    t - time index for meanV
%    ci - the 95% conf. interval around meanV
%    cispan -  the width of ci
%    prob - the probability of the model fit
%    QCspan - the 95th percentile minus the 50th percentile of the timestack
%             histogram, used to measure the amount of video "texture"
%    stdV - the width (std. dev.) of the energy in velocity spectrum
%    vAngle - orientation of the pixel array (radians)

%%  Criterios de Aceptación (Según Chickadel 2003)
NL=length(videoCurrentOut(200).InstVbar.meanI); %300 250
for j=1:si(1,2);
    if isstruct(videoCurrentOut(j).InstVbar)
        VQO=videoCurrentOut(j).InstVbar;   
        %NL=length(VQO.meanI);
        
        % (1) Greater than 90% significance of fit from the model skill
        CritAcep(j).CA(:,1)=[VQO.prob > P]';
        CritAcep(j).CountCA1=sum(CritAcep(j).CA(:,1));

        % (2) 95% confidence range of less than 0.2 m/s (+-0.1 m/s)
        CritAcep(j).CA(:,2)= [VQO.cispan < Lci]';
        CritAcep(j).CountCA2=sum(CritAcep(j).CA(:,2));

        % (3) Irange > 40
        CritAcep(j).CA(:,3)= [VQO.QCspan > Ir]';
        CritAcep(j).CountCA3=sum(CritAcep(j).CA(:,3));

        % Mean currents were considered valid if at least 10 of the 63 estimates
        % for each time series passed the above criteria (Aprox 16%)
        CritAcep(j).CountCA1(:,2)=(CritAcep(j).CountCA1(:,1)/NL);
        CritAcep(j).CountCA1(:,3)=(CritAcep(j).CountCA1(:,1)/NL)>0.16;

        CritAcep(j).CountCA2(:,2)=(CritAcep(j).CountCA2(:,1)/NL);
        CritAcep(j).CountCA2(:,3)=(CritAcep(j).CountCA2(:,1)/NL)>0.16;

        CritAcep(j).CountCA3(:,2)=(CritAcep(j).CountCA3(:,1)/NL);
        CritAcep(j).CountCA3(:,3)=(CritAcep(j).CountCA3(:,1)/NL)>0.16;

        Aux(1,1)=CritAcep(j).CountCA1(:,3);
        Aux(1,2)=CritAcep(j).CountCA2(:,3);
        Aux(1,3)=CritAcep(j).CountCA3(:,3);
        CritAcep(j).SeAcepta=sum(Aux)>2;
        
    else
        CritAcep(j).CA(:,1)=zeros(NL,1);
        CritAcep(j).CA(:,2)= zeros(NL,1);
        CritAcep(j).CA(:,3)= zeros(NL,1);
        
        CritAcep(j).CountCA1=0;
        CritAcep(j).CountCA1(:,2)=0;
        CritAcep(j).CountCA1(:,3)=0;
        
        CritAcep(j).CountCA2=0;
        CritAcep(j).CountCA2(:,2)=0;
        CritAcep(j).CountCA2(:,3)=0;
        
        CritAcep(j).CountCA3=0;
        CritAcep(j).CountCA3(:,2)=0;
        CritAcep(j).CountCA3(:,3)=0;
        
        CritAcep(j).SeAcepta=0;
    end

end

%% Output/Saving

% Output Name
oname=[NameV '_OutputChickadel_Tw' num2str(Tw) '_Ts' num2str(Ts)];

% OutPut Directory
odir='./Outputs/2-OutputChickadel';

% Save file
save([odir '/' oname ],'videoCurrentOut','CritAcep','errores')