function frameRect=makeRectSingleFramePracticum(I,xy,z, beta, lcp)
%Description:
%This practicum script rectifies a single image (I) using the specified real
%world coordinates (xy,z), the camera geometry or extrinsics (beta), and 
%the camera instrinsics (lcp).

%Inputs:
%I = image to be rectified
%xy = vector of xy information [xmin dx xmax ymin dy ymax]
%z = single elevation or grid of z elevations
%beta = geometry (extrinsics) from specified frame [camera X, camera Y, 
%camera Z, azimuth, tilt, roll] NOTE:all angles are in radians in your
%local coordinate system
%lcp = intrinsic values

%Outputs:
%frameRect = rectified image structure
%frameRect.x = x coordinates [1xN]
%frameRect.y = y coordinates [1xM]
%frameRect.I = rectified image [MxNx3]

%% organize indices
%I = double(I);
[NV,NU,NC] = size(I);
UI = [1:NU];
VI = [1:NV]';

%% define x,y,z grids
x = [xy(1):xy(2): xy(3)]; y = [xy(4):xy(5): xy(6)];
[X,Y] = meshgrid(x,y);

if length(z)==1
    xyz = [X(:) Y(:) repmat(z, size(X(:)))];
else
    xyz = [X(:) Y(:) z(:)];
end

%% Recall, Projection Matrix, P=KR[I|-C]
%Calculate P matrix
%define K matrix (intrinsics)
K = [lcp.fx 0 lcp.c0U;  
     0 -lcp.fy lcp.c0V;
     0  0 1];
%define rotation matrix, R (extrinsics)
R = angles2R(beta(4), beta(5), beta(6));
%define identity & camera center coordinates
IC = [eye(3) -beta(1:3)'];
%calculate P
P = K*R*IC;
%make P homogenous
P = P/P(3,4);   

%% Now, convert XYZ coordinates to UV coordinates
%convert xyz locations to uv coordinates
UV = P*[xyz'; ones(1,size(xyz,1))];

%homogenize UV coordinates (divide by 3 entry)
UV = UV./repmat(UV(3,:),3,1);

%convert undistorted uv coordinates to distorted coordinates
[Ud,Vd] = distort(UV(1,:),UV(2,:),lcp); 
UdVd = round([Ud; Vd]);%round to the nearest pixel locations
UdVd = reshape(UdVd,[],2); %reshape the data into something useable

%find the good pixel coordinates that are actually in the image (your grid
%might extend outside your image)
good = find(onScreen(UdVd(:,1),UdVd(:,2),NU,NV));


%% Finally, grab the RGB intensities at the indices you need and fill into your XYgrid
%preallocate final orthophoto
%preallocate final orthophoto
Irect = zeros(length(y),length(x),3);

% cycle through R,G,B intensities
for i = 1:NC    
    singleBandImage = I(:,:,i); %extract the frame
    rgbIntensities = interp2(UI, VI, double(singleBandImage),Ud,Vd); %extract the data at the pixel locations you need
    tempImage = Irect(:,:,i); %preallocate the orthophoto size based on your x,y 
    tempImage(good) = rgbIntensities(good); %fill your extracted pixels into the frame
    Irect(:,:,i) = tempImage; %put the frame into the orthophoto
end

frameRect.x=x;
frameRect.y=y;
frameRect.I=Irect;

%% Plot
imagesc(frameRect.x,frameRect.y,uint8(frameRect.I));
title('Rectified Image')
xlabel('X (m)');ylabel('Y (m)')
axis xy;axis image
end