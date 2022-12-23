% sample an example Aerielle movie.  Set this up as generic for other
% movies 

% The important things are:
%   - the init material below.  Stored as a structure so I can preserve
%       recordsfor each UAV analysis.
%   - a gcp file from a current survey
%   - 

clear
close all

addpath(genpath('.')) % adds folder dependencies 

 load('.\Outputs\OutputD-Part1.mat')
 
%% Modificación EstebanOV


cnt=1;
nStart=0;
h = waitbar(0,'Rectificando imágenes, espere porfavor...');
for clip = 1: NClips 
nStart=nStart+1;
    for i = nStart: Nf(clip) % eventually the length of fns
        cnt = cnt+1;


% read a new frame and find the new geometry, betas, and save in
% matrix.  Then sample both the stacks and build the image
% products, if needed.

I = imread([inputs.pnIn filesep clipFns(clip).name filesep fns{clip}(i).name]);

% Ejemplo:
% inputs.pnIn = 'C:\Codigos_Matlab\Support-Routines-Desarrollo1\StarterGuide\4-UAVArgus\demoMovies'
% clipFns(clip).name = 'demoClip1'
%fns{clip}(i).name = '00001001.jpg'

Ig = double(rgb2gray(I));


% findNewBeta(Ig,betas(cnt-1,~meta.globals.knownFlags), meta);

xyz = [meta.refPoints.xyz];
xyz = reshape(xyz(:),3,[])';
dUV = [meta.refPoints.dUV];
dUV = reshape(dUV(:),2,[])'; % Nose que sigbifica este dUV
thresh = [meta.refPoints.thresh]; %  Nose para que es esto (?????)

%%% [Ur,Vr, fail] = findCOMRefObj(Ig,xyz,beta,dUV,thresh,meta);
i;
ii=1;
fail = 0;               % flag to indicate no bright target found
minNGood = 4;           % fail if we don't find at least this # pixels
uv = round(findUVnDOF(betas(cnt-1,~meta.globals.knownFlags),xyz(ii,:),meta.globals));
URef = [uv(1)-dUV(1,1): uv(1)+dUV(1,1)];
VRef = [uv(2)-dUV(1,2): uv(2)+dUV(1,2)];
% you may not go off the edge!
VRef = VRef( find(VRef>0)); VRef = VRef( find(VRef<size(I,1)+1));
URef = URef( find(URef>0)); URef = URef( find(URef<size(I,2)+1));



I2 = I(VRef,URef);
[U,V] = meshgrid(URef,VRef);
good = find(I2>thresh(ii));
% if (length(good) < minNGood)
%     Ur = []; Vr = [];
%     fail = 1;
%     return
% end
Ur(ii) = mean(U(good));
Vr(ii) = mean(V(good));
% plot option
if meta.showFoundRefPoints
    figure(1+10);clf; colormap(gray)
    imagesc(URef,VRef,I2); axis image
    hold on
    plot(Ur(ii),Vr(ii),'r*')
    drawnow
end
du = round(Ur(ii) - uv(1));     % rough corrections to search guesses
dv = round(Vr(ii) - uv(2));



for ii = 2: size(xyz,1)
    uv = round(findUVnDOF(betas(cnt-1,~meta.globals.knownFlags),xyz(ii,:), meta.globals));
    uv = uv(:) + [du; dv];    
    URef = [uv(1)-dUV(ii,1): uv(1)+dUV(ii,1)];    %This is a problem if the search region heads out of view!!  (Negative indices!)
    VRef = [uv(2)-dUV(ii,2): uv(2)+dUV(ii,2)];  
    % you may not go off the edge! 
     URef = URef( find(URef>0)); URef = URef( find(URef<size(I,2)+1)); 
     VRef = VRef( find(VRef>0)); VRef = VRef( find(VRef<size(I,1)+1));
    I2 = I(VRef,URef);
    [U,V] = meshgrid(URef,VRef);
    Ur(ii) = mean(U(I2>thresh(ii)));
    Vr(ii) = mean(V(I2>thresh(ii)));
    % plot option
    if meta.showFoundRefPoints
        figure(ii+10);clf; colormap(gray)
        imagesc(URef,VRef,I2); axis image
        hold on
        plot(Ur(ii),Vr(ii),'r*')
        drawnow
    end
end
% if fail
%     beta6dof = [];
%     Ur = []; Vr = [];
%     return
% end



% Return findNewBeta(Ig,betas(cnt-1,~meta.globals.knownFlags), meta)

% find n-dof solution and expand to 6 dof final result
betaNew = nlinfit(xyz,[Ur(:); Vr(:)],'findUVnDOF',betas(cnt-1,~meta.globals.knownFlags));
beta6dof(find(meta.globals.knownFlags)) = meta.globals.knowns;
beta6dof(find(~meta.globals.knownFlags)) = betaNew;

% Para comprobar los Puntos de Referencia

if( meta.showInputImages == 1 )
    % show results in case debug is needed
    figure(2); clf; colormap(gray)
    imagesc(Ig); axis image
    hold on
    plot(Ur,Vr,'r*')
    uv = findUVnDOF(beta6dof, xyz, meta.globals);
    uv = reshape(uv,[],2);
    plot(uv(:,1),uv(:,2),'ko')
    hold off; drawnow;
end



beta1=beta6dof;
failFlag=fail;
%if failFlag         % deal with end of useable run
 %          break
  %   else
 betas(cnt,:) = beta1;
 
   %  end
% protect from running out of betas while images still exist
% with old meta file.
   %     if isnan(betas(cnt,1))
   %         break;
   %     end
   
   
   
   % 4. Perform the rectification
    % organize indices
    [NV,NU,NC] = size(I);
    Us = [1:NU];
    Vs = [1:NV]';

    % define x,y,z grids
        xy=inputs.rectxy;
        z=0;
        x = [xy(1):xy(2): xy(3)]; y = [xy(4):xy(5): xy(6)];
        [X,Y] = meshgrid(x,y);
        if length(z)==1
            xyz = [X(:) Y(:) repmat(z, size(X(:)))];
        else
            xyz = [X(:) Y(:) z(:)];
        end

    % Recall, Projection Matrix, P=KR[I|-C]
    %Calculate P matrix
    %define K matrix (intrinsics)
    K = [meta.globals.lcp.fx 0 meta.globals.lcp.c0U;  
         0 -meta.globals.lcp.fy meta.globals.lcp.c0V;
         0  0 1];
     

    %define rotation matrix, R (extrinsics)
    R = angles2R(beta1(4), beta1(5), beta1(6));
    %define identity & camera center coordinates
    IC = [eye(3) -beta1(1:3)'];
    %calculate P
    P = K*R*IC;
    %make P homogenous
    P = P/P(3,4);   
    

    % Now, convert XYZ coordinates to UV coordinates
    %convert xyz locations to uv coordinates
    UV = P*[xyz'; ones(1,size(xyz,1))];
    %homogenize UV coordinates (divide by 3 entry)
    UV = UV./repmat(UV(3,:),3,1);

    %convert undistorted uv coordinates to distorted coordinates
    [U,V] = distort(UV(1,:),UV(2,:),meta.globals.lcp); 
    UV = round([U; V]);%round to the nearest pixel locations
    UV = reshape(UV,[],2); %reshape the data into something useable

    %find the good pixel coordinates that are actually in the image 
    good = find(onScreen(UV(:,1),UV(:,2),NU,NV));
    %convert to indices
    ind = sub2ind([NV NU],UV(good,2),UV(good,1));

    % Finally, grab the RGB intensities at the indices you need and fill into your XYgrid
    %preallocate final orthophoto
    Irect = zeros(length(y),length(x),3);
    for iii = 1:NC    % cycle through R,G,B intensities
        singleBandImage = I(:,:,iii); %extract the frame
        rgbIndices = singleBandImage(ind); %extract the data at the pixel locations you need
        tempImage = Irect(:,:,iii); %preallocate the orthophoto size based on your x,y 
        tempImage(good) = rgbIndices; %fill your extracted pixels into the frame
        Irect(:,:,iii) = tempImage; %put the frame into the orthophoto
    end

    frameRect.x=x;
    frameRect.y=y;
    frameRect.I=Irect;
    Irect2=flipud(Irect);

    %Guardar Frames
    %Seleccionar carpeta que contiene carpeta Resultados
    imwrite(uint8(Irect2),'./Outputs/2-FramesRect/FrameRec_'+string(i)+'.jpg','jpg');
    save('./Outputs/2-FramesRect/FrameRec_'+string(i)+'.mat','Irect')
    
    %Guardar Metadata
    DataFrames.Coord.x(i,:)=x;
    DataFrames.Coord.y(i,:)=y;
    
    
    waitbar(i/Nf(clip),h)
    end
   
end
DataFrames.Time=dn;
waitbar(1,h,'El proceso ha finalizado correctamente')
clearvars -except DataFrames
save('./Outputs/DataFrames')

%
%   Copyright (C) 2017  Coastal Imaging Research Network
%                       and Oregon State University

%    This program is free software: you can redistribute it and/or  
%    modify it under the terms of the GNU General Public License as 
%    published by the Free Software Foundation, version 3 of the 
%    License.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see
%                                <http://www.gnu.org/licenses/>.

% CIRN: https://coastal-imaging-research-network.github.io/
% CIL:  http://cil-www.coas.oregonstate.edu
%
%key UAVProcessingToolbox
%

