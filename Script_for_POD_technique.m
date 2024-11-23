% Script for extracting dominant wave characteristics using POD technique. 
% This technique is featured in the paper:

% Bheeroo, V., Bae, S. B., Lee, M.- J., Socolofsky, S. A., & Chang, K.- A. 
% (2024). "Using unmanned aerial systems for observations of water wave 
% characteristics". Experiments in Fluids. 

% This script is used to extract the dominant wave characteristics for the
% wavefield presented in Sec. 3.1 of the paper. 
tic
clear
clc
close all; 
%% Extracting frames from video 
vidReader = VideoReader('UAS-Video_Gal_entrance_channel_20240215.MOV');
jj = 1;
frameIndex = 1;

while hasFrame(vidReader)
    frames_r = readFrame(vidReader);

    if mod(frameIndex, 3) == 1  % Save every 3rd frame so that the original video
                                % is downsampled from 24 fps to 8 fps
        im{jj,1} = frames_r(:,:,3); %Extracting blue channel
        fprintf('%d\n', jj);
        jj = jj + 1;
    end

    frameIndex = frameIndex + 1;
end
conv_px = 0.0273; %Image resolution in meters 
dt = 1/8 ;        %Sampling rate in seconds
%% Choosing a crop from the field of view for computational efficiency 
rect = load('cropped_region.mat');  %Region of crop used in paper
rect = rect.rect; 
sample_image = im{1, 1}; 
cropped_sample = imcrop(sample_image, rect);
[height, width, num_channels] = size(cropped_sample);  % Dimensions of cropped image
im_crop = zeros(length(im), height, width); %Pre-allocating cropped image matrix
for jj=1:length(im) 
    im_crop(jj,:,:) = imcrop(im{jj,1},rect); 
    fprintf('%d\n',jj)
end
Nx = size(im_crop,3); %Pixel dimensions in x 
Ny = size(im_crop,2); %Pixel dimensions in y
Nt = size(im_crop,1); %Number of image snapshots
x = linspace(0,Nx*conv_px,Nx); %x array
y = linspace(0,Ny*conv_px,Ny); %y array
t = linspace(0,Nt*dt,Nt); %t array
%% Applying POD to images and processing mode and mode coefficient
clc
[U,S,V] = pod(single(im_crop)); %See paper for definitions of U,S,V

mode_2 = squeeze(U(2,:,:));   %Mode 2 represents dominant wave signal
mode_2_coefficient = V(:,2); 

mode_2_autocorrelation = xcorr2(mode_2); %2D-Autocorrelation of mode 2
lags_x = linspace((-Nx-1)*conv_px,(Nx-1)*conv_px,size(mode_2_autocorrelation,2)); 
lags_y = linspace((-Ny-1)*conv_px,(Ny-1)*conv_px,size(mode_2_autocorrelation,1)); 
mode_2_autocorrelation_x = mode_2_autocorrelation(floor(size(mode_2_autocorrelation,1)/2),:);
mode_2_autocorrelation_y = mode_2_autocorrelation(:,floor(size(mode_2_autocorrelation,2)/2)); 
[mode_2_coefficient_autocorrelation,lags_t] = xcorr(mode_2_coefficient,'normalized'); %1D Autocorrelation of mode 2 coefficients
lags_t = lags_t*dt; 

[pks_x,lcs_x] = findpeaks(mode_2_autocorrelation_x,lags_x,'SortStr','descend'); 
[pks_y,lcs_y] = findpeaks(mode_2_autocorrelation_y,lags_y,'SortStr','descend'); 
[pks_t,lcs_t] = findpeaks(mode_2_coefficient_autocorrelation,lags_t,'SortStr','descend'); 

L_x = lcs_x(2); 
L_y = lcs_y(2); 
Wave_period = abs(lcs_t(2))
Wavelength = (L_x*L_y)/sqrt(L_x^2 + L_y^2)
Wave_direction = atan2(L_x,L_y) * (180/pi) %Wave direction from positive y-axis

toc