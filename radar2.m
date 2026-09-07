%% 三维 Sample × Chirp × Frame FMCW 仿真（exp 解析信号版本）
clear; close all; clc;

%% 雷达系统参数设置
fc= 77e9;             % 雷达工作频率
c = 3e8;              % 光速

%% FMCW波形参数设置
maxR = 200;           % 最大探测距离
rangeRes = 1;         % 距离分辨率
maxV = 20;            % 最大速度（未显式使用）
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
Mix = zeros(Nr, Nd, Nf);   % 复数混频信号（已为基带）

%% 逐帧生成数据（解析信号，exp 形式）
for f = 1:Nf
    % 本帧总时长 Nd*Tchirp，对应 Nr*Nd 个采样点
    t_frame = linspace(0, Nd*Tchirp, Nr*Nd); 
    
    % 帧起始距离（目标按恒速匀移）
    r0_frame = r0 + v0*(f-1)*(Nd*Tchirp); 
    
    % 目标往返时延
    r_t = r0_frame + v0*t_frame;     % 距离随时间
    td  = 2*r_t/c;                    % 延迟
    
    % 发射/接收相位
    phi_tx = 2*pi*( fc*t_frame + 0.5*slope.*t_frame.^2 );
    phi_rx = 2*pi*( fc*(t_frame-td) + 0.5*slope.*(t_frame-td).^2 );
    
    % 解析信号（复指数）
    Tx = exp(1j*phi_tx);
    Rx = exp(1j*phi_rx);
    
    % 复数共轭解调，直接得到拍频基带
    mix_frame = Tx .* conj(Rx);    % 已是基带复信号
    
    % 转成 Nr × Nd 并写入本帧
    Mix(:,:,f) = reshape(mix_frame, Nr, Nd);
end

%% 一些坐标轴
doppler_axis = linspace(0,128,Nd)*vres;
range_axis   = linspace(0,256,Nr)*rangeRes;

%% 1️⃣ 第一张图：第一帧的第一个 chirp 的 距离FFT 曲线
range_fft_firstchirp = abs( fft( Mix(:,1,1) ) );
figure;
plot(range_fft_firstchirp);
xlabel('距离（频率bin）');
ylabel('幅度');
title('第一个chirp的FFT结果 (第一帧)');

%% 2️⃣ 第二张图：第一帧所有chirp 的 距离维FFT谱矩阵（Nr × Nd）
sig_fft1 = fft( Mix(:,:,1), [], 1 );   % 对第1维（采样/距离）做FFT
sig_fft  = abs(sig_fft1);

figure;
mesh(sig_fft);
ylabel('距离（频率bin）');
xlabel('chirp脉冲数');
zlabel('幅度');
title('距离维FFT结果 (第一帧)');

%% 3️⃣ 第三张图：第一帧的 距离-多普勒谱（先距FFT，再沿chirp做多普勒FFT）
sig_fft2 = fft( sig_fft1, [], 2 );     % 对第2维（chirp）做FFT
RDM = abs(sig_fft2);

figure;
mesh(doppler_axis, range_axis, RDM);
xlabel('多普勒通道'); ylabel('距离通道'); zlabel('幅度');
title('速度维FFT 距离-多普勒谱 (第一帧)');

%% 4️⃣ 新增图：第一个 chirp 的 距离×时间（帧）热力图
range_time_map = zeros(Nr, Nf);
for f = 1:Nf
    range_time_map(:,f) = abs( fft( Mix(:,1,f) ) );
end

figure;
mesh(1:Nf, range_axis, range_time_map);
set(gca,'YDir','normal');
xlabel('时间（帧）');
ylabel('距离');
title('第一个chirp的 距离×时间 图');
colorbar;
