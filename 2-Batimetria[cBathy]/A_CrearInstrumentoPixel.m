%% Creación Instrumento Pixel Grid

% Section 1: User Input:  Saving Information
clear
res=5; %Resolucion [m]
oname=['15V3_' num2str(res) 'm'];

%  Enter the directory where the instrument file will be saved.
odir='./Outputs/1-PixelInstruments/';

%% Section 2: User Input:  Collection Information

%Imagenes Rectificadas
imageDirectory{1}= 'C:\Memoria2020\Resultados\Rectificados\15V3\MAT';

%Info Imagenes Rectificadas

DireccionDF='./Inputs/DataFrames_15V3';

% DirGuard = uigetdir([],'Seleccione carpeta donde se ubican los frames rectificados'); %Resultados
% clipFns = dir(DirGuard);

%% Section 3: User Input: Manual Entry of Time and Elevation

t={};
zFixedCam={};

%% Section 4: User Input: Instrument Information


localFlag=1;

%  Grid
%  Enter the following parameters for a grid. Note, dx and dy do not
%  need to be equal.
pixInst(1).type='Grid';
pixInst(1). dx =res;
pixInst(1). dy =res;
pixInst(1). xlim =[80 400];
pixInst(1). ylim =[-300 300];
pixInst(1).z={}; % Leave empty if you would like it interpolated from input
% Z grid or zFixedCam. If entered here it is assumed constant
% across domain and in time.

%% Section 5: Load Files 

load(DireccionDF)
x=DataFrames.Coord.x(1,:);
y=DataFrames.Coord.y(1,:);
sx=length(x); %size x
sy=length(y); %size y
X=repmat(x,sy,1);
Y=repmat(y,sx,1)';
Z=zeros(sy,sx);

% Load List of Collection Images
L{1}=string(ls(imageDirectory{1}));


%% Section 6: Initialize Pixel Instruments

% Make XYZ Grids/Vectors from Specifications.
[pixInst]=pixInstPrepXYZ(pixInst);

% Assign Z Elevations Depending on provided parameters.
% If pixInst.z is left empty, assign it correct elevation
for p=1:length(pixInst)
    if isempty(pixInst(p).z)==1
        if isempty(zFixedCam)==1 % If a Time-varying z is not specified
            pixInst(p).Z=interp2(X,Y,Z,pixInst(p).X,pixInst(p).Y); % Interpret from given grid
        end
        
        if isempty(zFixedCam)==0 % If a time varying z is specified
            pixInst(p).Z=pixInst(p).X.*0+zFixedCam(1); % Assign First Value for First Image, spatially constant elevation
        end
    end
end

%% Section 7:  Loop for Collecting Pixel Instrument Data.
tic

%Arreglo de indices
XY(:,1) = X(:); XY(:,2) = Y(:); %Coordenadas imagen original
XYins(:,1) = pixInst.X(:); XYins(:,2) = pixInst.Y(:); %Coordenadas Instrumento Grid
ism = ismember(XY,XYins,'rows');
ind = find(ism == 1);

 for j=1:(length(L{1}(:))-2)

    % Load Imagen Rectificada
    Irgb=load(strcat(imageDirectory{1}, '\FrameRec_',num2str(j),'.mat'));
    %I=I1.Irect;

    %  Loop for Each Pixel Instrument
    for p=1:length(pixInst)
        
        % Check if a time varying Z was specified. If not, wil just use constant Z
        % specified or interpolated from grid in Section 4 and 6 respectively.
        if isempty(pixInst(p).z)==1
            if isempty(zFixedCam)==0 % If a time varying z is specified
                pixInst(p).Z=pixInst(p).X.*0+zFixedCam(j); % Assign First Value for First Image, spatially constant elevation
            end
        end

        
        % Pull RGB Pixel Intensities
        %[Irgb]= I;
        
        
        % Convert To Grayscale
        
        [Igray]=rgb2gray(uint8(Irgb.Irect));
        clear Irgb
        Ig=reshape(Igray(ind),size(pixInst.X));
        
        
        % If First frame, initialize pixInst structure entries
        if j==1
            pixInst(p).Igray=Ig;
%            pixInst(p).Irgb=Irgb;
        end
        
        % If not First frame, tack on as last dimension (time).
        if j~=1
            s=size(Ig);
            nDim= length(find(s~=1)); % Finds number of actual dimension is =1 (transects) or 2 (grid)nDim
            
            % For Gray Scale it is straight forward
            pixInst(p).Igray=cat(nDim+1,pixInst(p).Igray,Ig); % Add on last dimension (third if Grid, second if transects)
            
            % For RGB it is trickier since MATLAB always likes rgb values in
            % third dimension.
            % If a GridGrid Add in the fourth dimension
            if nDim==2
%(Por Memoria)                pixInst(p).Irgb=cat(nDim+2,pixInst(p).Irgb,Irgb); % Add on last dimension (Always Fourth)
            end
            % If a Transect Grid Add in the second dimension
            if nDim==1
%(Por Memoria)                pixInst(p).Irgb=cat(2,pixInst(p).Irgb,Irgb);
            end
            
        end
   end
    

    
end
toc
%% Section 9: Save Instruments and Meta Data

% Save metaData and Grid Data

rectMeta.imageDirectory=imageDirectory;
rectMeta.imageNames=L;

save(strcat(odir, '/',oname, '_pixInst'), 'pixInst','rectMeta','t')



