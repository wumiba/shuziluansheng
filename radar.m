%% 三维 Sample × Chirp × Frame FMCW 仿真
clear; close all; clc;

%% 雷达系统参数设置
fc= 77e9;             % 雷达工作频率
c = 3e8;              % 光速

%% FMCW波形参数设置
maxR = 200;           % 最大探测距离
rangeRes = 1;         % 距离分辨率
maxV = 20;            % 最大速度
B = c/(2*rangeRes);   % 带宽
Nd=128;               % chirp数
Nr=256;               % 采样点数
Nf=32;                % 帧数（新增）
Tchirp = 5 * 2 * maxR/c;  
endle_time=2.5e-6;    

vres=(c/fc)/(2*Nd*(Tchirp+endle_time)); % 速度分辨率

%% 用户目标参数
r0 = rangeRes*90;     % 初始距离
v0 = vres*5;          % 速度

%% 调频斜率
slope = B / Tchirp; 

%% 数据存储 (三维矩阵 Nr × Nd × Nf)
Mix = zeros(Nr, Nd, Nf);

%% 逐帧生成数据
for f = 1:Nf
    t_frame = linspace(0, Nd*Tchirp, Nr*Nd); 
    Tx = zeros(1, length(t_frame));
    Rx = zeros(1, length(t_frame));
    mix_frame = zeros(1, length(t_frame));
    
    r0_frame = r0 + v0*(f-1)*(Nd*Tchirp); % 每帧起始距离
    
    for i = 1:length(t_frame)
        r_t = r0_frame + v0*t_frame(i);
        td  = 2*r_t/c;
        
        Tx(i) = cos(2*pi*(fc*t_frame(i) + (slope*t_frame(i)^2)/2));
        Rx(i) = cos(2*pi*(fc*(t_frame(i)-td) + (slope*(t_frame(i)-td)^2)/2));
        mix_frame(i) = Tx(i).*Rx(i);
    end
    
    Mix(:,:,f) = reshape(mix_frame, Nr, Nd);
end

%% 1️⃣ 第一个chirp的距离FFT
MixIQ = zeros(Nr,Nd);
for i=1:Nd
    MixIQ(:,i) = hilbert(Mix(:,i,1)); % 希尔伯特变换
end

sig_fft1 = zeros(Nr,Nd);
for k=1:Nd
    sig_fft1(:,k)=fft(MixIQ(:,k));
end
sig_fft = abs(sig_fft1);

figure;
plot(sig_fft(:,1));
xlabel('距离（频率）');
ylabel('幅度')
title('第一个chirp的FFT结果 (第一帧)')

%% 2️⃣ 距离FFT结果谱矩阵（第一帧）
figure;
mesh(sig_fft);
ylabel('距离（频率）');
xlabel('chirp脉冲数')
zlabel('幅度')
title('距离维FFT结果 (第一帧)')

%% 3️⃣ 速度维FFT（第一帧）
sig_fft2 = zeros(Nr,Nd);
for k=1:Nr
    sig_fft2(k,:)=fft(sig_fft1(k,:));
end

RDM = abs(sig_fft2);
doppler_axis = linspace(0,128,Nd)*vres;
range_axis = linspace(0,256,Nr)*rangeRes;

figure;
mesh(doppler_axis,range_axis,RDM);
xlabel('多普勒通道'); ylabel('距离通道'); zlabel('幅度（dB）');
title('速度维FFT 距离-多普勒谱 (第一帧)');

%% 4️⃣ 新增图：第一个chirp的距离FFT，显示 距离 × 时间（帧）
range_time_map = zeros(Nr,Nf);
for f = 1:Nf
    % 取每帧的第一个chirp
    chirp_data = hilbert(Mix(:,1,f)); 
    range_fft = abs(fft(chirp_data));
    range_time_map(:,f) = range_fft;
end

figure;
mesh(1:Nf, range_axis, range_time_map);
xlabel('时间（帧）');
ylabel('距离');
title('第一个chirp的 距离×时间 图');
colorbar;
