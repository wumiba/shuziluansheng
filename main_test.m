%{
clear;
a = load('D:\FMCWdata\3D非视距代码\NLOS_3D_60_1\adcData64.mat');
b = load('D:\FMCWdata\3D非视距代码\NLOS_3D_60_2\adcData64.mat');


%}
if ~((exist('a', 'var') ~= 0) || (exist('b', 'var') ~= 0) || (exist('c', 'var') ~= 0))
    a = load('D:\FMCWdata\3D非视距代码\NLOS_3D_60_2\adcData64.mat');
end
    

clearvars -except a CDMheat CDMbreath;
close;
%a = load('D:\FMCWdata\3D非视距代码\NLOS_3D_60_2\adcData64.mat');
adcData64 = a.adcData64;
para = f_parameter64();
adcData64 = reshape(adcData64, para.AntNum, para.adcsamples,para.chirploops, para.datalength);

%adcData64 = adcData64 - repmat(mean(adcData64,4),[1 1 1 1200]);

chirp64 = squeeze(adcData64(1, :, :, :));

% sumchirp = squeeze(chirp64(:,1,:));

sumchirp = squeeze(sum(chirp64,2));



%{
c = load('E:\非视距雷达生命体征数据\NLOS_Data\NLOS_4D_45_1\NLOS_4D_45_1.mat');
c = load('E:\非视距雷达生命体征数据\NLOS_64G_matlab_data\RX1.mat');
variableNames = fieldnames(c);
channel_data = variableNames{1};
adcData64 = permute(c.(channel_data), [3 4 1 2 5]);
chirp64 = squeeze(adcData64(1, 1 , :, :, 9:end));
para = f_parameter64();
para.datalength = 1190;
para.chirploops = 4;
sumchirp = chirp64;
sumchirp = squeeze(sum(chirp64,2));
%}

%测试多普勒测速显示
%{
    figure;
    FFTSize = para.adcsamples;
    fRangeresol = para.rangemax/FFTSize;
    Rvec = (0:FFTSize-1)*fRangeresol;    

    x = (-para.chirploops/2 : para.chirploops/2 - 1)*(para.vmax/para.chirploops);
    y = Rvec;
    RDdata = squeeze(adcData64(1 , 1 , : , : , 3));
    mesh( x , y , abs(fft2(RDdata))) ;
    xlabel("速度 (m/s)");
    ylabel("距离 (m)");

    axis([-0.3 0.3 0 10]);

%}

%Rangefft去均值后测速
%{
    figure;
    FFTSize = para.adcsamples;
    fRangeresol = para.rangemax/FFTSize;
    Rvec = (0:FFTSize-1)*fRangeresol;   
    RDdata = squeeze(adcData64(1 , 1 , : , : , :));
    D_RangeFFT = zeros(FFTSize,para.chirploops,para.datalength);
    D_RangeFFT = fft(RDdata , FFTSize , 1);
    D_mean_data = mean(D_RangeFFT,3);
    D_RangeFFTcor = D_RangeFFT - D_mean_data;
    D_RangeFFTcor = squeeze(D_RangeFFTcor(: , : , 5));


    x = (-para.chirploops/2 : para.chirploops/2 - 1)*(para.vmax/para.chirploops);
    y = Rvec;

    mesh( x , y , fftshift(abs(fft(D_RangeFFTcor,4,2)) , 2));
    xlabel("速度 (m/s)");
    ylabel("距离 (m)");

    axis([-0.3 0.3 0 10]);

%}

%测试角度速度
%{
    figure;
    FFTSize = para.adcsamples;
    fRangeresol = para.rangemax/FFTSize;
    Rvec = (0:FFTSize-1)*fRangeresol;

    test_fftsize = 1024;
    x = (-para.chirploops/2 : para.chirploops/2 - 1)*(para.vmax/para.chirploops);
    y = asind((0 : test_fftsize - 1)/test_fftsize*2-1);
    RDdata = squeeze(adcData64( : , 1 , : , 1));
    mesh( x , y , fftshift(abs(fft2(RDdata , 1024 , 64))));
    xlabel("速度");
    ylabel("角度");

%}


%测试角度距离
%{
    
    figure;
    FFTSize = para.adcsamples;
    fRangeresol = para.rangemax/FFTSize;
    Rvec = (0:FFTSize-1)*fRangeresol;

    test_fftsize = 1024;
    NewfRangeresol = fRangeresol / (test_fftsize / para.adcsamples);
    x = ((0 : test_fftsize - 1)) * NewfRangeresol;
    y = asind((0 : test_fftsize - 1)/test_fftsize*2-1);
    RDdata = squeeze(sum(adcData64(: , : , : , : , 1),4));
    mesh( x , y , fftshift(abs(fft2(RDdata , test_fftsize , test_fftsize)) , 1));
    xlabel("距离");
    ylabel("角度");
    axis([6 8 -100 100 0 130000]);

%}


%{

非相干积累
chirp64 = squeeze(sum(abs(adcData64).^2,1));

sumchirp = squeeze(sum(abs(chirp64).^2,2));

%}

% sumchirp = diff(sumchirp , 1 , 2);




%% Obtain range profile
FFTSize = para.adcsamples;
fRangeresol = para.rangemax/FFTSize;
Rvec = (0:FFTSize-1)*fRangeresol;


RangeFFT = zeros(FFTSize,para.datalength);
for ii =1:para.datalength
    RangeFFT(:,ii) = fft(hanning(length(sumchirp(:,ii))).*sumchirp(:,ii),FFTSize);
end

% for tt = 1:270
%     RangeFFT(tt,:) = RangeFFT(120,:);
% end
% 
% for pp = 420:FFTSize
%     RangeFFT(pp,:) = RangeFFT(120,:);
% end

%进行真正的慢时间维度非相干积累
%{

FFTSize = para.adcsamples;
fRangeresol = para.rangemax/FFTSize;
Rvec = (0:FFTSize-1)*fRangeresol;
RangeFFT = zeros(FFTSize,para.datalength);
for ii =1:para.datalength
    for jj = 1 : 64
        to_fft = squeeze(chirp64(:,jj,:));
        to_abs = fft(hanning(length(to_fft(:,ii))).*to_fft(:,ii),FFTSize);
        RangeFFT(:,ii) = RangeFFT(:,ii) + abs(to_abs).^2;
    end
end
%}

% FFTSize = para.adcsamples;
% fRangeresol = para.rangemax/FFTSize;
% Rvec = (0:FFTSize-1)*fRangeresol;
% RangeFFT = zeros(FFTSize,para.datalength);
% for ii =1:para.datalength
%     for jj = 1 : 4
%         to_fft = squeeze(chirp64(:,jj,:));
%         to_abs = fft(hanning(length(to_fft(:,ii))).*to_fft(:,ii),FFTSize);
%         RangeFFT(:,ii) = RangeFFT(:,ii) + abs(to_abs).^2;
%     end
% end

%% DC offset correction
mean_data = mean(RangeFFT,2);
RangeFFTcor = RangeFFT - mean_data;

figure;
mesh(Rvec,(0:para.datalength-1)*50*1e-3,abs(RangeFFT)'/max(abs(RangeFFT(:))))
view([0,90])
colorbar;
xlabel("距离 (m)");
ylabel("时间 (s)");
title("反正切相位序列结果图");
xlim([0 7]);
%axis([6 8 1 (para.datalength-1)*50*1e-3]);

% [position, threshold] = CFAR(abs(RangeFFTcor(:,1)'), 0.15, 10, 3);
% figure;
% plot(Rvec,threshold);
% hold on;grid on;
% plot(Rvec,abs(RangeFFTcor(:,1)'));




%% （原有代码保持不变，直到RangeFFTcor计算结束）
%% MTI处理（一次对消）
% 沿慢时间维度（datalength）做相邻帧差分
RangeFFT_MTI = diff(RangeFFT, 1, 2);  % 1阶差分，沿第2维度（datalength）

%% 显示MTI后的结果
figure;
% 调整时间轴（差分后帧数减少1）
time_axis_MTI = (0:size(RangeFFT_MTI,2)-1) * 50e-3;  


mesh(Rvec, time_axis_MTI, abs(RangeFFT_MTI)' / max(abs(RangeFFT_MTI(:))));
view([0, 90]);
colorbar;
xlabel("距离 (m)");
ylabel("慢时间 (s)");
title("MTI处理后的距离-时间图");
axis([0 (FFTSize-1)*fRangeresol 0 (size(RangeFFT_MTI,2)-1)*50e-3]);




















MaxIndex = 100;
MaxIndex2 = 93;

It = zeros(1,para.datalength);
Qt = zeros(1,para.datalength);

It2 = zeros(1,para.datalength);
Qt2 = zeros(1,para.datalength);

for ii = 1: para.datalength
    It(ii) = real(RangeFFT(MaxIndex,ii));
    Qt(ii) = imag(RangeFFT(MaxIndex,ii));
    
    It2(ii) = real(RangeFFTcor(MaxIndex2,ii));
    Qt2(ii) = imag(RangeFFTcor(MaxIndex2,ii));
end

arctan_index2=atan2(It2,Qt2);
plot((0:para.datalength-1)*50*1e-3,arctan_index2);
xlabel('时间 (s)');
ylabel('相位 (rad)');
title('正切相位序列时域图');

N = para.datalength;
for n=1:N
    phi1(n)=0;
    phi2(n)=0;
    for k=2:n
        phi1(n)=phi1(n)+(It(k)*(Qt(k)-Qt(k-1))-Qt(k)*(It(k)-It(k-1)))/(It(k).^2+Qt(k).^2);
        phi2(n)=phi2(n)+(It2(k)*(Qt2(k)-Qt2(k-1))-Qt2(k)*(It2(k)-It2(k-1)))/(It2(k).^2+Qt2(k).^2);
    end
end




angle=diff(phi2);


    for ii=3:para.datalength-1
        angleDenoised(ii)=f_phaseDenoise(angle(ii-2),angle(ii-1),angle(ii),0.3);
    end

%angle=improvedRemoveImpulseNoise(angleDenoised);
%end


%{
%plot phase data and frequency spectrum
figure;
subplot(121)
plot((0:length(angleDenoised)-1)*50e-3,angleDenoised,'LineWidth',2)
xlabel("时间 (s)");
ylabel("相位 (rads)");
title("差分相位序列时序图");

subplot(122)
plot((-length(angleDenoised)/2:length(angleDenoised)/2-1)*20/length(angleDenoised),abs(fftshift(fft(angleDenoised)))/max(abs(fftshift(fft(angleDenoised)))),'LineWidth',2)
axis([0 10 0 1.1]);
xlabel("频率 (Hz)");
ylabel("幅度");
title("差分相位序列频率图");
%}

sc = 3e8/(4*pi*7.8995e+10);
fs =para.fps; %Sampling frequency

%
%% Signal separation through Variational Mode Decomposition变分模态分解
k = 7;  %Number of modes
alpha = 1000;    %penalty factor
tol = 1e-6; %Absolute tolerance
[imf,IMFs_fft,~] = VMD(angleDenoised,alpha,0,k,0,2,tol);

IMFs_fft=abs(IMFs_fft)';
M=size(IMFs_fft,2);
freqs = (-M/2:M/2-1)*(fs/M);


%{
figure;
for k = 1:size(imf,1)
subplot(size(imf,1),1,k)
plot((0:1000-1)*(60/1000),imf(k,:))
title(['IMF',num2str(k)])
xlabel('时间 (s)');
ylabel('幅度');
end



figure;
for k = 1:size(imf,1)
subplot(size(imf,1),1,k)
plot(freqs,IMFs_fft(k,:));
axis([0 10 0 max(abs(IMFs_fft(k,:)))])
title(['IMF',num2str(k),' 频谱图'])
xlabel('频率 (Hz)');
ylabel('幅度');
end


figure;
for k = 1:4
subplot(4,1,k)
plot((0:1000-1)*(60/1000),imf(k,:))
title(['IMF',num2str(k)])
xlabel('时间 (s)');
ylabel('幅度');
end



figure;
for k = 1:4
subplot(4,1,k)
plot(freqs,IMFs_fft(k,:));
axis([0 10 0 max(abs(IMFs_fft(k,:)))])
title(['IMF',num2str(k),' 频谱图'])
xlabel('频率 (Hz)');
ylabel('幅度');
end

figure;
for k = 5:7
subplot(4,1,k-4)
plot((0:1000-1)*(60/1000),imf(k,:))
title(['IMF',num2str(k)])
xlabel('时间 (s)');
ylabel('幅度');
end



figure;
for k = 5:7
subplot(4,1,k-4)
plot(freqs,IMFs_fft(k,:));
axis([0 10 0 max(abs(IMFs_fft(k,:)))])
title(['IMF',num2str(k),' 频谱图'])
xlabel('频率 (Hz)');
ylabel('幅度');
end
%}


fshift = freqs;
resp_range = [0.1 0.5];  % Respiration frequency range
hb_range = [0.8 2];      % Heart beat frequency range

resp_idx = find(fshift >= resp_range(1) & fshift <= resp_range(2));
hb_idx = find(fshift >= hb_range(1) & fshift <= hb_range(2));

resp_amplitudes = zeros(1, size(imf, 1));  % Preallocate array for respiration amplitudes
hb_amplitudes = zeros(1, size(imf, 1));    % Preallocate array for heart beat amplitudes

for i = 1:size(imf, 1)
    hb_amplitudes(i) = max(IMFs_fft(i, hb_idx));
end
[~, hb_idx] = max(hb_amplitudes);

for i = 1:size(imf, 1)
    resp_amplitudes(i) = max(IMFs_fft(i, resp_idx));
end

[~, resp_idx] = max(resp_amplitudes);
resp_IMF = imf(resp_idx, :);
breath_fre = IMFs_fft(resp_idx,:);
hb_IMF = imf(hb_idx, :);
heart_fre = IMFs_fft(hb_idx,:);

%disp(['Respiratory IMf: imf',num2str(resp_idx),' Heart beat IMf: imf',num2str(hb_idx)])
% spectral analysis for VMD
%
% Find breathing frequency
[pks, locs] = findpeaks(breath_fre);
[~, idx] = max(pks);
breath_count = (fs * ((para.datalength-1)/2-(locs(idx)-1))/(para.datalength-1)) * 60;

% Find heart beat frequency
% Compute power spectral density of heart signal
[Pxx, f] = pwelch(hb_IMF, [], [], [], fs);

% Keep only frequencies in range [0.8, 2] Hz
mask = f >= 0.8 & f <= 2;
Pxx_masked = Pxx(mask);
f_masked = f(mask);

% Find maximum in PSD
[~, idx] = max(Pxx_masked);

% Convert peak location to Hz and compute heart rate
heart_count = f_masked(idx) * 60;


disp(['BR: ',num2str(breath_count),' bpm','    HR: ',num2str(heart_count),' bpm'])



%% Breathing signals filtering呼吸信号滤波
%{
COE1=Butterworth_IIR; %filterdesigner generated filter
save coe1.mat COE1;
breath_data = filter(COE1,angleDenoised); 
%breath_datam1 = filter(COE1,angleDenoised(2,:)); 
%breath_datap1 = filter(COE1,angleDenoised(3,:)); 

%{
figure;
plot((0:length(breath_data)-1)*50e-3,breath_data*sc,'LineWidth',2);
xlabel('时间 (s)');
ylabel('位移 (m)');
title('呼吸信号时域波形');
%}


%% Heart signal filtering心跳信号滤波
COE2=Butterworth_IIR2;
save coe2.mat COE2;
heart_data = filter(COE2,angleDenoised);
%heart_datam1 = filter(COE2,angleDenoised(2,:));
%heart_datap1 = filter(COE2,angleDenoised(3,:));

%{
figure;
plot((0:length(heart_data)-1)*50e-3,heart_data*sc,'LineWidth',2);
xlabel('时间(s)');
ylabel('位移 (m)');
title('心跳信号时域波形');
%}

%% Spectral estimation谱估计
N1=length(breath_data);
fshift = (-N1/2:N1/2-1)*(fs/N1); % zero-centered frequency scale
breath_fre = abs(fftshift(fft(breath_data)));              %--FFT
%}
%{
figure;
plot(fshift,breath_fre,'LineWidth',2);
xlabel('frequency (Hz)');
ylabel('Amplitude');
title('Breathing signal FFT');
axis([0 10 0 max(breath_fre)+1])


% Find breathing frequency
[pks, locs] = findpeaks(breath_fre);
[~, idx] = max(pks);
breath_count = (fs * ((para.datalength-1)/2-(locs(idx)-1))/(para.datalength-1)) * 60; %Breathing rate

N1=length(heart_data);
fshift = (-N1/2:N1/2-1)*(fs/N1); % zero-centered frequency
heart_fre = abs(fftshift(fft(heart_data))); 
%}
%{
figure;
plot(fshift,heart_fre,'LineWidth',2);
xlabel('frequency (Hz)');
ylabel('Amplitude');
title('Heart beat FFT');
axis([0 10 0 max(heart_fre)+1])
%}
% Heart rate estimation
%{
[pks, locs] = findpeaks(heart_fre);
[~, idx] = max(pks);
heart_count = (fs * ((para.datalength-1)/2-(locs(idx)-1))/(para.datalength-1)) * 60;;%Heart rate estimation

disp(['BR: ',num2str(breath_count),' bpm','    HR: ',num2str(heart_count),' bpm'])

%}
%}




%%  Perform CEEMDAN on phase signal对相位信号执行CEEMDAN
%{
Nstd = 2*var(angleDenoised); % noise standard deviation
NR = 500;   %number of realisations
MaxIter = 7; %Max number of sifting operations
SNRFlag = 2;    %Signal to Noise Ratio flag
[modes,its] = ceemdan(angleDenoised,Nstd, NR, MaxIter, SNRFlag);

M = length(modes);

freqs = (-M/2:M/2-1)*(fs/M);
IMFs_fft = zeros(size(modes));
for i=1:size(modes,1)
    IMFs_fft(i,:) = abs(fftshift(fft(modes(i,:))));
end

% plot the IMFs and their frequency spectra
%{

figure;
for i=1:size(modes,1)
    subplot(size(modes,1),2,2*i-1);
    plot(modes(i,:));
    title(['IMF ' num2str(i)]);
    xlabel('时间 (s)');
    ylabel('幅度');
    subplot(size(modes,1),2,2*i);
    plot(freqs, IMFs_fft(i,:));
    title(['IMF ' num2str(i) ' 频谱']);
    xlabel('频率 (Hz)');
    ylabel('幅度');
    axis([0 10 0 max(IMFs_fft(i,:))])
end

figure;
for i=1:size(modes,1)/2
    subplot(size(modes,1)/2,2,2*i-1);
    plot((0:length(angleDenoised)-1)*50e-3,modes(i,:));
    title(['IMF ' num2str(i)]);
    xlabel('时间 (s)');
    ylabel('幅度');
    subplot(size(modes,1)/2,2,2*i);
    plot(freqs, IMFs_fft(i,:));
    title(['IMF ' num2str(i) ' 频谱']);
    xlabel('频率 (Hz)');
    ylabel('幅度');
    axis([0 10 0 max(IMFs_fft(i,:))])
end
figure;
for i=size(modes,1)/2:size(modes,1)-1
    subplot(size(modes,1)/2,2,2*(i-4)-1);
    plot((0:length(angleDenoised)-1)*50e-3,modes(i,:));
    title(['IMF ' num2str(i)]);
    xlabel('时间 (s)');
    ylabel('幅度');
    subplot(size(modes,1)/2,2,2*(i-4));
    plot(freqs, IMFs_fft(i,:));
    title(['IMF ' num2str(i) ' 频谱']);
    xlabel('频率 (Hz)');
    ylabel('幅度');
    axis([0 10 0 max(IMFs_fft(i,:))])
end
figure;
i=10;
subplot(4,2,1);
    plot((0:length(angleDenoised)-1)*50e-3,modes(i,:));
    title(['IMF ' num2str(i)]);
    xlabel('时间 (s)');
    ylabel('幅度');
    subplot(4,2,2);
    plot(freqs, IMFs_fft(i,:));
    title(['IMF ' num2str(i) ' 频谱']);
    xlabel('频率 (Hz)');
    ylabel('幅度');
    axis([0 10 0 max(IMFs_fft(i,:))])


%}
%}





%数字滤波呼吸心跳
%{

diff_phase = angleDenoised;
fs = 20; % 采样频率（Hz）

T = length(diff_phase) / fs; % 数据总时长（秒）

% 设计低频带通滤波器 (0.1~0.5 Hz)
low_freq_band = [0.1, 0.5]; 
[b_low, a_low] = butter(5, low_freq_band / (fs/2), 'bandpass'); % 4阶Butterworth
filtered_low = filtfilt(b_low, a_low, diff_phase);

% 设计高频带通滤波器 (0.8~2 Hz)
high_freq_band = [0.8, 2];
[b_high, a_high] = butter(4, high_freq_band / (fs/2), 'bandpass');
filtered_high = filtfilt(b_high, a_high, diff_phase);

% 绘制时域结果
t = (0:length(diff_phase)-1) / fs; % 时间轴（秒）
figure;
plot(t, diff_phase);
title('原始差分相位序列');
xlabel('时间 (s)');

figure;
plot(t, filtered_low);
title('呼吸信号时域图');
xlabel('时间 (s)');
ylabel("幅度");

figure;
plot(t, filtered_high);
title('心跳信号时域图');
xlabel('时间 (s)');
ylabel("幅度");

figure;
plot((-length(filtered_low)/2:length(filtered_low)/2-1)*20/length(filtered_low),abs(fftshift(fft(filtered_low)))/max(abs(fftshift(fft(filtered_low)))),'LineWidth',2)
axis([0 5 0 1.1]);
xlabel("频率 (Hz)");
ylabel("幅度");
title("呼吸信号频谱图");

figure;
plot((-length(filtered_high)/2:length(filtered_high)/2-1)*20/length(filtered_high),abs(fftshift(fft(filtered_high)))/max(abs(fftshift(fft(filtered_high)))),'LineWidth',2)
axis([0 5 0 1.1]);
xlabel("频率 (Hz)");
ylabel("幅度");
title("心跳信号频谱图");

[~ ,num_ofheart] = findSignalPeaks(hb_IMF,20,'on')


huitu.data1 = abs(fftshift(fft(filtered_high)))/max(abs(fftshift(fft(filtered_high))));
huitu.data1_x = (-length(filtered_low)/2:length(filtered_low)/2-1)*20/length(filtered_low);
huitu.name1 = '数字滤波';
huitu.data2 = abs(((heart_fre)))/max(abs(((heart_fre))));
huitu.data2_x = freqs;
huitu.name2 = 'CEEMDAN';
huitu.x_label = '频率 (Hz)';
huitu.y_label = '幅度';
huitu.data_title = '心跳信号频谱对比图';
jiahuatu = huatu(huitu);
xlim([0 5]);


huitu.data1 = abs(fftshift(fft(filtered_low)))/max(abs(fftshift(fft(filtered_low))));
huitu.data1_x = (-length(filtered_low)/2:length(filtered_low)/2-1)*20/length(filtered_low);
huitu.name1 = '数字滤波';
huitu.data2 = abs(((breath_fre)))/max(abs(((breath_fre))));
huitu.data2_x = freqs;
huitu.name2 = 'CEEMDAN';
huitu.x_label = '频率 (Hz)';
huitu.y_label = '幅度';
huitu.data_title = '呼吸信号频谱对比图';
jiahuatu = huatu(huitu);
xlim([0 5]);




%CDM与VDM对比
huitu.data1 = abs(((CDMheat)))/max(abs(((CDMheat))));
huitu.data1_x = (-length(filtered_low)/2:length(filtered_low)/2-1)*20/length(filtered_low);
huitu.name1 = 'CEEMDAN';
huitu.data2 = abs(((heart_fre)))/max(abs(((heart_fre))));
huitu.data2_x = freqs;
huitu.name2 = 'VDM';
huitu.x_label = '频率 (Hz)';
huitu.y_label = '幅度';
huitu.data_title = '心跳信号频谱对比图';
jiahuatu = huatu(huitu);
xlim([0 5]);


huitu.data1 = abs(((CDMbreath)))/max(abs(((CDMbreath))));
huitu.data1_x = (-length(filtered_low)/2:length(filtered_low)/2-1)*20/length(filtered_low);
huitu.name1 = 'CEEMDAN';
huitu.data2 = abs(((breath_fre)))/max(abs(((breath_fre))));
huitu.data2_x = freqs;
huitu.name2 = 'VDM';
huitu.x_label = '频率 (Hz)';
huitu.y_label = '幅度';
huitu.data_title = '呼吸信号频谱对比图';
jiahuatu = huatu(huitu);
xlim([0 5]);



%VDM结果图
figure;
plot((0:1000-1)*(60/1000),hb_IMF)
title(['心跳信号时域图'])
xlabel('时间 (s)');
ylabel('幅度');


figure;
plot(freqs,heart_fre);
xlim([0 10])
title(['心跳信号频谱图'])
xlabel('频率 (Hz)');
ylabel('幅度');

figure;
plot((0:1000-1)*(60/1000),resp_IMF)
title(['呼吸信号时域图'])
xlabel('时间 (s)');
ylabel('幅度');


figure;
plot(freqs,breath_fre);
xlim([0 10])
title(['呼吸信号频谱图'])
xlabel('频率 (Hz)');
ylabel('幅度');













%CDM结果图
figure;
plot((0:1199-1)*(60/1200),hb_IMF)
title(['心跳信号时域图'])
xlabel('时间 (s)');
ylabel('幅度');


figure;
plot(freqs,heart_fre);
xlim([0 10])
title(['心跳信号频谱图'])
xlabel('频率 (Hz)');
ylabel('幅度');

figure;
plot((0:1199-1)*(60/1200),resp_IMF)
title(['呼吸信号时域图'])
xlabel('时间 (s)');
ylabel('幅度');


figure;
plot(freqs,breath_fre);
xlim([0 10])
title(['呼吸信号频谱图'])
xlabel('频率 (Hz)');
ylabel('幅度');

%}

