% Background/Foreground Separation via Dynamic Mode Decomposition (DMD)
% Johnathon R Barhydt
%
% This code takes a video and uses the DMD modes with the smallest
% time-variance to determine background, from which a foreground detection 
% filter is made and masked with the original video.
%
% specs for this run:
% two videos from cell phone, about 10 seconds long. first video is looking
% out the window of a coffee shop, the camera is still. second video is
% held by hand, and shows dogs running around a dog park. background is
% more shaky in the video.
%
clear all; close all; clc
% video name
vid_name='DCIM_0001.mp4';
% rescale video
a=1/8;
% import video
v = VideoReader(vid_name);
n_frames = 1;
while hasFrame(v)
    vid = readFrame(v);
    video(:,:,:,n_frames) = imresize(vid, [a*v.Height a*v.Width]);
    n_frames = n_frames+1;
end
% video dimensions
im_h=size(video,1);
im_w=size(video,2);
% remove fencepost frame
n_frames=n_frames-1;
% save video as matrix
save('video.mat','video');

%%
% import saved video (redundant if already performed above)
vid_mat='video.mat';
% get camera data
video = importdata(vid_mat);
%%
% create colummated data matrix, one snapshot per column
for i=1:size(video,4)
    image = im2double( imresize( rgb2gray(video(:,:,:,i)) ,[im_h im_w] ) );
    X(:,i)= reshape(image,[im_h*im_w,1]);
end

% send data matrix to my DMD function
vid_filter = mat2fg(X);
%% show fg video filter
for i=1:n_frames
    frame=reshape( vid_filter(:,i), [im_h im_w]);
    imshow(mat2gray(frame)); drawnow
end
%% show fg filtered video
% triplicate filter for r,g,b channels
A=cat(3,vid_filter,vid_filter,vid_filter);
A=permute(A,[1 3 2]);
A=reshape( A,[im_h im_w 3 n_frames]);

% filtered video, with BG removed!
f_vid = double(video(:,:,:,1:end-1)).*double(A);

for i=1:n_frames
    imshow(uint8(f_vid(:,:,:,i))); drawnow
end

%% For FUN! this plots ALL oscillation modes and their evolution in space
% figure(1)
% set(gcf, 'Position', get(0, 'Screensize'));
% set(gca,'Ylim',[-.00001 .00001],'Xlim',[-.00001 .00001],'Color','k')
% for i=1:449   
%     plot(b_e_wt(i,:),'w.-')
%     %plot(exp(w(i)*i),'ko')
%     set(gca,'Ylim',[-.00001 .00001],'Xlim',[-.00001 .00001],'Color','k')
%     pause(0.08)
%     drawnow
% end




