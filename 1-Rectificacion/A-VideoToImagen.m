% Video a Imagenes


clear
addpath(genpath('.'))



%% SECCION 1: Inputs 


nVid='DJI_0001.mov'; %Nombre Archivo Video 
fpsR=2; % Razón de extracción de frames, ¿cuántos frames deseo extraer por segundo(fps)?

%Others
SaveImages = 1;  % Do you want to save individual frames from the video? 1=yes 0=no. 
% CollectFrames = 0; % Do you want to collect frames into array and save? 1=yes 0=no. 

%Direcion Output
DirOut='Outputs\1-Frames\';

%% SECCION 2: Parametros


% Construct a VideoReader object associated with the file  
[pathstr,name,ext]=fileparts(nVid);  
obj = VideoReader(nVid); 
 
% obj contains information about the video, e.g.:
 vidHeight = obj.Height; %height of image in pixels
 vidWidth = obj.Width; %width of image in pixels
 vidFrameRate = obj.FrameRate;  % Maybe you should round this to the nearest int.
 
% Resample at desired framerate:
 numberOfFrames = floor(obj.Duration / (1/obj.FrameRate)); %Much faster calculation from given values.
 framesToRead = 1: round(round(vidFrameRate)/fpsR) :numberOfFrames; %%Esta parte modifiqué, le agregué un round
 
% if CollectFrames ==1 % Only do this if you intend to save the resampled frames as an array
%     allFrames = zeros(vidHeight, vidWidth, 3, length(framesToRead));  % Preallocate array to store frames.
% end
   
   
%% SECCION 3: MAIN LOOP

tic

 h = waitbar(0,'Espere un momento mientras se extraen los frames...');
 for k = 1:length(framesToRead)
     
  %if mod(k,10)==0; disp(k); end % counter to indicate something happening. 
   frameIdx = framesToRead(k);
   currentFrame = read(obj,frameIdx);


   if SaveImages == 1   % Save the frames of interest as jpg images
     combinedString = sprintf([name '_' '%05d.jpg'],k);
     imwrite(currentFrame,[DirOut filesep combinedString]);
   end
 
%    if CollectFrames ==1  %Gather the frames of interest into a giant 4D array (RGB and time)
%      if k==1
%        % cast the all frames matrix appropriately
%        allFrames = cast(allFrames, class(currentFrame)); % This step takes long time!
%      end
%      allFrames(:,:,:,k) = currentFrame;
%    end

   waitbar(k/length(framesToRead),h)
 end
 
 
% if CollectFrames ==1  % SAVE THE ARRAY.
%     % Save the array of frames for future use:
%     save([pathstr filesep name '_resampled_at_' num2str(fpsR) 'fps.mat'], 'allFrames')
% end
   
   waitbar(1,h,'El proceso ha finalizado correctamente')
 
toc

%% Output/Saving

DataVideo.Duration=obj.Duration;

DataVideo.Path=obj.Path;
DataVideo.BitsPerPixel=obj.BitsPerPixel;
DataVideo.FrameRate=obj.FrameRate;
DataVideo.Height=obj.Height;
DataVideo.NumFrames=obj.NumFrames;
DataVideo.VideoFormat=obj.VideoFormat;
DataVideo.Width=obj.Width;

VideoName=obj.Name;
VideoName=VideoName(1:end-4)
DataVideo.Name=VideoName;
% Save file

save(['./Outputs/DataVideo-' VideoName],'DataVideo')
save('./Outputs/DataAnalisis','DataVideo')
 
 
 
 