%% platform IWR1642EVM+DCA1000
%% vital sign monitoring
%% single target: one male adult, sitting 1.4m from radar, for 51.2 seconds (1024frame)
%% please check /experimental_settings.jpg
%% ========================================================================
clc;
clear all;
close all;
%% =========================================================================
%% chirp settings
%% please check /transceiver_settings.png
numADCSamples = 200; % number of ADC samples per chirp
numADCBits = 16;     % number of ADC bits per sample
numRX = 4;           % number of receivers
numTX = 2;           % number of transceiver
numLanes = 2;        % do not change. number of lanes is always 2
isReal = 0;          % set to 1 if real only data, 0 if complex data0
chirpLoop = 2;       % number of chirps per frame
Fs=4e6;              % sampling rate of ADC
c=3*1e8;             % speed of light
ts=numADCSamples/Fs; % ADC sampling time 
slope=70e12;         % slope of the chirp 
B_valid =ts*slope;   % bandwith of the chirp
detaR=c/(2*B_valid); % range resolution

%% read data file
Filename = 'normal_breath.bin';  % normal_breath.bin or fast_breath.bin
fid = fopen(Filename,'r');
adcDataRow = fread(fid, 'int16');
if numADCBits ~= 16
    l_max = 2^(numADCBits-1)-1;
    adcDataRow(adcDataRow > l_max) = adcDataRow(adcDataRow > l_max) - 2^numADCBits;
end
fclose(fid);

fileSize = size(adcDataRow, 1);
PRTnum = fix(fileSize/(numADCSamples*numRX));
fileSize = PRTnum * numADCSamples*numRX;
adcData = adcDataRow(1:fileSize);
% real data reshape, filesize = numADCSamples*numChirps
if isReal
    numChirps = fileSize/numADCSamples/numRX;
    LVDS = zeros(1, fileSize);
    %create column for each chirp
    LVDS = reshape(adcData, numADCSamples*numRX, numChirps);
    %each row is data from one chirp
    LVDS = LVDS.';
else
    numChirps = fileSize/2/numADCSamples/numRX;     %complex data include real and imaginary part
    LVDS = zeros(1, fileSize/2);
    %combine real and imaginary part into complex data
    %read in file: 2I is followed by 2Q
    counter = 1;
    for i=1:4:fileSize-1
        LVDS(1,counter) = adcData(i) + sqrt(-1)*adcData(i+2);
        LVDS(1,counter+1) = adcData(i+1)+sqrt(-1)*adcData(i+3); counter = counter + 2;
    end
    % create column for each chirp
    LVDS = reshape(LVDS, numADCSamples*numRX, numChirps);
    %each row is data from one chirp
    LVDS = LVDS.';
end

%% reshape data
adcData = zeros(numRX,numChirps*numADCSamples);
for row = 1:numRX
    for i = 1: numChirps
        adcData(row, (i-1)*numADCSamples+1:i*numADCSamples) = LVDS(i, (row-1)*numADCSamples+1:row*numADCSamples);
    end
end

retVal= reshape(adcData(1, :), numADCSamples, numChirps); %select the data from receivers1,the data of each chirp is stored in a row

process_adc=zeros(numADCSamples,numChirps/2);

for nchirp = 1:2:numChirps  %select the data of transceiver1
    process_adc(:, (nchirp-1)/2+1) = retVal(:,nchirp);
end

	
%% range FFT of the first chirp (For display purposes only)
figure;
rangefft_chirp_1 = db(abs(fft(process_adc(:,1))));
% DC interference mitigation
for i = 1:3 
     rangefft_chirp_1(i) = rangefft_chirp_1(i) * 0.8;
end
plot((1:numADCSamples)*detaR,rangefft_chirp_1,'LineWidth',1.5);
xlabel('range£¨m£©');
ylabel('Amp(dB)');
title('rangeFFT');

%% range FFT
RangFFT = 256;
fft_data_last = zeros(1,RangFFT); 
range_max = 0;
adcdata = process_adc;
numChirps = size(adcdata, 2);

for nchirp = 1:2:numChirps 
   adcdata(:, (nchirp-1)/2+1) = retVal(:,nchirp);  % there are 2 chirps in a frame, choose the first chirp
end
numChirps = numChirps/2;


fft_data = fft(adcdata,RangFFT); 
fft_data = fft_data.';

for ii=1:numChirps-1                % Sliding cancellation 
     fft_data(ii,:) = fft_data(ii+1,:)-fft_data(ii,:);
end

fft_data_abs = abs(fft_data);

fft_data_abs(:,1:10)=0; % remove DC component

real_data = real(fft_data);
imag_data = imag(fft_data);


% extract phase of every range bin(frequency), using arctangent
for i = 1:numChirps
    for j = 1:RangFFT  
        angle_fft(i,j) = atan2(imag_data(i, j),real_data(i, j));
    end
end

% Range-bin selecting 
% a simplest algorithm: find the range bin(frequency) with the largest total amp, to be the range bin of the human target
% this Range-bin selecting algorithm is very rough and sensitive to clutter, such as furniture or obstacles.
% the Range-bin is correctly selected only when no strong clutter existed in the monitoring area.
% Range-bin selecting algorithm is one of the problems that I'm going to focus on in PhD study.
for j = 1:RangFFT
    if((j*detaR)<1.7 &&(j*detaR)>1.2)      % this algorithm only works in clear area
        for i = 1:numChirps  % sum up amp of each range bin(frequency)
            fft_data_last(j) = fft_data_last(j) + fft_data_abs(i,j);
        end
        
        if ( fft_data_last(j) > range_max)
            range_max = fft_data_last(j);
            max_num = j;  
        end
    end
end 

%% extract phase from selected range bin
angle_fft_last = angle_fft(:,max_num);

%% phase unwrapping.
%% Since the phase value is between [-pi, pi], and we need phase unwrapping to obtain the actual displacement curve,
%% whenever the phase difference between continuous values is greater than pi or less than -pi, 
%% the phase unwrapping is obtained by subtracting or adding 2pi from the phase.

n = 1;
for i = 1+1:numChirps
    diff = angle_fft_last(i) - angle_fft_last(i-1);
    if diff > pi
        angle_fft_last(i:end) = angle_fft_last(i:end) - 2*pi;
        n = n + 1;
    elseif diff < -pi
        angle_fft_last(i:end) = angle_fft_last(i:end) + 2*pi;  
    end
end

%% phase difference. This will help to enhance the heartbeat signal and eliminate the existing phase drift
angle_fft_last2=zeros(1,numChirps);
for i = 1:numChirps-1
    angle_fft_last2(i) = angle_fft_last(i+1) - angle_fft_last(i);
    angle_fft_last2(numChirps)=angle_fft_last(numChirps)-angle_fft_last(numChirps-1);
end 

figure;
plot(angle_fft_last2);
xlabel('time(frame)');
ylabel('phase');
title('phase waveform');

%% Bandpass Filter 0.1-0.6hz for respiration signal
fs =20; % sampling rate of vital sign signal (rate of the frame)

COE1=respiration_filter;
breath_data = filter(COE1,angle_fft_last2); 

figure;
plot(breath_data,'LineWidth',1.5);
xlabel('time(frame)');
ylabel('amp');
title('respiration waveform in time domain');

%% fft to obtain the spectrum of the respiration signal
N1=length(breath_data);
fshift = (-N1/2:N1/2-1)*(fs/N1); % zero-centered frequency
breath_fre = abs(fftshift(fft(breath_data)));             



breath_fre_max = 0; % respiration rate estimation
for i = 1:length(breath_fre) 
    if (breath_fre(i) > breath_fre_max)    
        breath_fre_max = breath_fre(i);
        breath_index=i;
    end
end

breath_count =(fs*(numChirps/2-(breath_index-1))/numChirps)*60; % respiration per minute

%% Bandpass Filter 0.8-2hz, for heartbeat signal
%% heartbeat signal filtered by bandpass filter, include 2th or 3th harmonics interference of respiration signal
%% the true heartbeat signal is easily submergered by harmonics when respiration movement is strong and fast.
%% more sophisticated methods to mitigate respiration harmonics in heartbeat signal will be explored in the future study
COE2=heartbeat_filter;
heart_data = filter(COE2,angle_fft_last2); 
figure;
plot(heart_data,'LineWidth',1.5);
xlabel('time(frame)');
ylabel('amp');
title('heartbeat waveform in time domain');

N1=length(heart_data);
fshift = (-N1/2:N1/2-1)*(fs/N1); % zero-centered frequency
heart_fre = abs(fftshift(fft(heart_data))); 


heart_fre_max = 0; 
for i = 1:length(heart_fre)/2 
    if (heart_fre(i) > heart_fre_max)    
        heart_fre_max = heart_fre(i);
        if(heart_fre_max<1e-2) % 
            heart_index=1025;
        else
            heart_index=i;
        end
    end
end
heart_count =(fs*(numChirps/2-(heart_index-1))/numChirps)*60;% heartbeat per minute

% 1024 frames£¬51.2s£¬

disp(['respiration per minute',num2str(breath_count),'  heartbeat per minute',num2str(heart_count)])

%% END &thank YOU !
