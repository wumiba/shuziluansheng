% combined_chest_heart_sim.m
% 将胸腔 (breath) 与心脏 (heartbeat) 两个模型合并在一份脚本中
% 每个模型保留自己原始的参数、形变公式与 STL 文件
% 最终在同一场景中绘图并计算 FMCW 雷达回波（两模型贡献叠加）

clear; close all; clc;

%% ---------------------- STL 文件名（按需修改） ----------------------
chest_stl = 'd.stl';   % 胸腔 STL（若你原来用 d.stl 可改回）
heart_stl = '心脏2_5.stl';   % 心脏 STL（若你原来用 b.stl 可改回）

%% ==================== === 胸腔（来源：你给出的“胸腔”代码） ====================
% 我将胸腔部分变量前缀为 chest_ 以和 heart 区分，保留原始参数
% ---------------------- 呼吸用户参数（胸腔） ----------------------
T_in       = 0.8;      % 吸气上升阶段时长 [s]
T_plateau  = 0.4;      % 平台期（吸气或呼气的维持）[s]
T_out      = 0.8;      % 呼气下降阶段时长 [s]
T_b        = T_in + T_plateau + T_out + T_plateau;   % 完整呼吸周期
T_h_chest  = T_b;        % 用呼吸周期覆盖原心跳周期（其余代码仍用 T_h）
fps        = 25;         % 帧率
tau_h_chest    = T_h_chest / 3;
sigma_h_chest  = T_h_chest / 5;
d_h_chest      = 0.002;       % 最大形变位移量 [m]
direction_chest = [-0.5, -1.0, 0.0];

% 心脏影响范围（这是胸腔脚本里的“heart region”）
x_min_chest = -0.13; x_max_chest = 0.13;
y_min_chest = -0.2;  y_max_chest = 0;
z_min_chest = -0.2;  z_max_chest = 0.2;

daxiao = 1; % 缩放系数（胸腔原来有）
x_min_chest = x_min_chest * daxiao;
x_max_chest = x_max_chest * daxiao;
y_min_chest = y_min_chest * daxiao;
y_max_chest = y_max_chest * daxiao;
z_min_chest = z_min_chest * daxiao;
z_max_chest = z_max_chest * daxiao;
d_h_chest = d_h_chest * daxiao;

% 动画长度（共用 fps）
total_seconds = 8;
total_frames = round(total_seconds * fps);

% chest stl 路径（使用 chest_stl 变量）

use_vertex_normals_chest = true; % chest 原脚本 true

%% ==================== === 心脏（来源：你给出的“心脏”代码） ====================
% 前缀 heart_
T_h_heart    = 0.8;        % 心跳周期（秒）
fps_heart    = 25;         % 帧率（与胸腔保持一致）
tau_h_heart  = T_h_heart / 3;
sigma_h_heart = T_h_heart / 5;
d_h_heart     = 0.003;     % 最大形变位移量 [m]
direction_heart = [-0.5, -1.0, 0.0];

% 心脏影响范围（心脏脚本里的）
x_min_heart = -0.016; x_max_heart = 0.03;
y_min_heart = -0.02;  y_max_heart = 0.01;
z_min_heart = -0.022; z_max_heart = 0.0075;


daxiao_heart = 2.5;
x_min_heart = x_min_heart*daxiao_heart; x_max_heart = x_max_heart*daxiao_heart;
y_min_heart = y_min_heart*daxiao_heart;  y_max_heart = y_max_heart*daxiao_heart;
z_min_heart = z_min_heart*daxiao_heart; z_max_heart = z_max_heart*daxiao_heart;
d_h_heart  = d_h_heart*daxiao_heart;

% 动画保持 total_seconds, total_frames 相同
% stl 文件 heart_stl

use_vertex_normals_heart = true; % 心脏原脚本 false (沿固定方向扩张)

%% ==================== === 雷达参数（默认采用胸腔代码的雷达参数） ===
% 注意：必须选择一组用于仿真的雷达参数。为了保持与原来胸腔脚本一致，
% 我把胸腔里的雷达参数设置为主参数（不修改其值）。
% fc = 77e9;             % 雷达工作频率 [Hz]
% c = 3e8;               % 光速 [m/s]
% 
% Nchirp = 2;              % chirp数（每个FMCW帧内的chirp数）
% Nsample = 64;           % 每个chirp的采样点数
% maxR = 20;            % 最大探测距离 [m]
% rangeRes = maxR/Nsample;          % 距离分辨率 [m]
% maxV = 20;             % 最大速度 [m/s]
% 
% B = c / (2 * rangeRes); % 带宽 [Hz]
% Nframe = total_frames;   % 帧数（与动画帧对应）
% 
% Tchirp = 2.5 * 2 * maxR / c;  % chirp持续时间 [s]
% idle_time = 5e-6;    % idle时间 [s]
% vres = (c / fc) / (2 * Nchirp * (Tchirp + idle_time)); % 速度分辨率 [m/s]
% 
% slope = B / Tchirp;    % 调频斜率 [Hz/s]

fc = 77.00e9;             % 雷达工作频率 [Hz]
c = 3e8;                  % 光速 [m/s]
B = 3.99034e9;            % 带宽 [Hz]
Nchirp = 2;              % 每帧 chirp 数
Nsample = 64;            % 每个 chirp 的采样点数
slope = 70.006e12;        % 调频斜率 [Hz/s]
Tchirp = B / slope;       % chirp 持续时间 [s]，由 B 和 slope 计算
sample_rate = 8000e3;     %采样速率
chirp_time = Nsample/sample_rate;  %有效chirptime
idle_time = 5e-6;         % chirp 间空闲时间 [s]
Nframe = total_frames; % 总帧数（原 datalength）
maxR = (sample_rate*c)/(2*slope);    % 最大探测距离 [m]
% rangeRes = maxR/Nsample; % 距离分辨率 [m]
value_B = chirp_time*slope;
rangeRes = c/(2*value_B)
lambda = c / fc;
dx = lambda / 2;  % Rx 之间的间隔：λ/2
N_rx = 8;


% 雷达位置（以 chest 脚本为主），你可改为任何值
p_radar = [0, -5, 0]; % 雷达位置 [m]
theta_0 = 0;           % 方位角参考
phi_0 = 0;             % 俯仰角参考

epsilon = 1e-12;       % 数值稳定性

% Mix 存储（Nsample × Nchirp × Nframe）
Mix = zeros(Nsample, Nchirp, Nframe);

%% ==================== GPU 初始化（全局） ====================
if gpuDeviceCount > 0
    gpu = gpuDevice(1);
    disp('✅ 使用GPU加速: ' + string(gpu.Name));
else
    error('❌ 未检测到兼容GPU。请检查硬件和MATLAB Parallel Computing Toolbox。');
end

%% ==================== 读取两个 STL 并分别预处理 ====================
% 使用原始 read_stl_generic 函数（在脚本末尾）
[Vertices_chest, Faces_chest] = read_stl_generic(chest_stl);
[Vertices_heart,  Faces_heart ] = read_stl_generic(heart_stl);

Vertices_chest = double(Vertices_chest);
Faces_chest = double(Faces_chest);
Vertices_heart = double(Vertices_heart);
Faces_heart = double(Faces_heart);

% 保持原始顶点索引（但合并场景时需要偏移）
numVerts_chest = size(Vertices_chest,1);
numVerts_heart = size(Vertices_heart,1);

% 计算各自顶点法向量（与原代码一致）
% chest normals
numFacesC = size(Faces_chest,1);
n_sum_c = zeros(numVerts_chest,3);
for f = 1:numFacesC
    i1 = Faces_chest(f,1); i2 = Faces_chest(f,2); i3 = Faces_chest(f,3);
    v1 = Vertices_chest(i1,:); v2 = Vertices_chest(i2,:); v3 = Vertices_chest(i3,:);
    nf = cross(v2 - v1, v3 - v1);
    nf = nf / (norm(nf) + epsilon);
    n_sum_c(i1,:) = n_sum_c(i1,:) + nf;
    n_sum_c(i2,:) = n_sum_c(i2,:) + nf;
    n_sum_c(i3,:) = n_sum_c(i3,:) + nf;
end
vertex_normals_chest = n_sum_c ./ (sqrt(sum(n_sum_c.^2,2)) + epsilon); % Nx3

% heart normals
numFacesH = size(Faces_heart,1);
n_sum_h = zeros(numVerts_heart,3);
for f = 1:numFacesH
    i1 = Faces_heart(f,1); i2 = Faces_heart(f,2); i3 = Faces_heart(f,3);
    v1 = Vertices_heart(i1,:); v2 = Vertices_heart(i2,:); v3 = Vertices_heart(i3,:);
    nf = cross(v2 - v1, v3 - v1);
    nf = nf / (norm(nf) + epsilon);
    n_sum_h(i1,:) = n_sum_h(i1,:) + nf;
    n_sum_h(i2,:) = n_sum_h(i2,:) + nf;
    n_sum_h(i3,:) = n_sum_h(i3,:) + nf;
end
vertex_normals_heart = n_sum_h ./ (sqrt(sum(n_sum_h.^2,2)) + epsilon); % Nx3

% 将法向量放到GPU（后续可能被两边共享）
vertex_normals_chest_gpu = gpuArray(vertex_normals_chest);
vertex_normals_heart_gpu = gpuArray(vertex_normals_heart);



%% ==================== 分别计算形变幅度 A_i（胸腔与心脏） ====================
% chest A (保留胸腔原来计算方法，注意 chest 原代码把 sqrt_fz = fz not sqrt)
e_const = exp(1);
A_chest = zeros(numVerts_chest,1);
for i=1:numVerts_chest
    x = Vertices_chest(i,1); y = Vertices_chest(i,2); z = Vertices_chest(i,3);
    if ~(x >= x_min_chest && x <= x_max_chest && y >= y_min_chest && y <= y_max_chest && z >= z_min_chest && z <= z_max_chest)
        continue;
    end
    fz = - ( (z - z_min_chest) * (z - z_max_chest) ) / ((z_max_chest - z_min_chest)^2);
    fz = max(fz, 0);
    fy = - ( (y - y_min_chest) * (y - y_max_chest) ) / ((y_max_chest - y_min_chest)^2);
    fy = max(fy, 0);
    sqrt_fz = fz; % 注意：胸腔原代码这里没有 sqrt
    sqrt_fy = sqrt(fy);
    ratio = (x - x_min_chest) / (x_max_chest - x_min_chest);
    ratio = min(max(ratio,0),1);
    log_term = log( e_const + (1 - e_const) * ratio );
    A_chest(i) = d_h_chest * sqrt_fz * 1 * log_term;
end
A_chest_gpu = gpuArray(A_chest);

% heart A (保持心脏原来的计算)
A_heart = zeros(numVerts_heart,1);
for i=1:numVerts_heart
    x = Vertices_heart(i,1); y = Vertices_heart(i,2); z = Vertices_heart(i,3);
    if ~(x >= x_min_heart && x <= x_max_heart && y >= y_min_heart && y <= y_max_heart && z >= z_min_heart && z <= z_max_heart)
        continue;
    end
    fz = - ( (z - z_min_heart) * (z - z_max_heart) ) / ((z_max_heart - z_min_heart)^2);
    fz = max(fz, 0);
    fy = - ( (y - y_min_heart) * (y - y_max_heart) ) / ((y_max_heart - y_min_heart)^2);
    fy = max(fy, 0);
    sqrt_fz = sqrt(fz);
    sqrt_fy = sqrt(fy);
    ratio = (x - x_min_heart) / (x_max_heart - x_min_heart);
    ratio = min(max(ratio,0),1);
    log_term = log( e_const + (1 - e_const) * ratio );
    A_heart(i) = d_h_heart * sqrt_fz * sqrt_fy * log_term;
end
A_heart_gpu = gpuArray(A_heart);

%% ==================== 确定形变方向（分别） ====================
% chest delta_max
if use_vertex_normals_chest
    delta_max_chest_gpu = A_chest_gpu .* vertex_normals_chest_gpu;
else
    dir_c = direction_chest / norm(direction_chest);
    dir_c_gpu = gpuArray(dir_c);
    delta_max_chest_gpu = A_chest_gpu .* dir_c_gpu;
end

% heart delta_max
if use_vertex_normals_heart
    delta_max_heart_gpu = A_heart_gpu .* vertex_normals_heart_gpu;
else
    dir_h = direction_heart / norm(direction_heart);
    dir_h_gpu = gpuArray(dir_h);
    delta_max_heart_gpu = A_heart_gpu .* dir_h_gpu;
end

%% ==================== 合并顶点与面（用于绘图与统一雷达处理） ====================
% 合并顶点
Vertices_all = [Vertices_chest; Vertices_heart];
Faces_all = [Faces_chest; Faces_heart + size(Vertices_chest,1)];
vertex_normals_all = [vertex_normals_chest; vertex_normals_heart];

% 合并 delta_max（顶点数相应拼接）
delta_max_all_gpu = [delta_max_chest_gpu; delta_max_heart_gpu];

% 将总体顶点推到 GPU，初始化 V_prev
Vertices_all_gpu = gpuArray(Vertices_all);
V_prev_gpu = Vertices_all_gpu;
p_radar_gpu = gpuArray(p_radar);

% 顶点法向量 GPU（用于 Lambertian）
vertex_normals_all_gpu = gpuArray(vertex_normals_all);

%% ==================== FMCW 时间轴（GPU） ====================
t_frame_fmcw_gpu = gpuArray(linspace(0, Nchirp * Tchirp, Nsample * Nchirp));  % [s]

%% ==================== 绘图初始化 ====================
fig = figure('Name','Chest + Heart Animation with Radar Echo','Color',[1 1 1]);
ax = axes('Parent', fig);
p = patch('Vertices', Vertices_all, 'Faces', Faces_all, ...
    'FaceColor',[0.8 0.3 0.3], 'EdgeColor','none', ...
    'FaceLighting','gouraud', 'AmbientStrength',0.3, ...
    'SpecularStrength',0.2, 'Parent', ax);
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
material shiny; camlight headlight; view(3); grid on;
vc = mean(Vertices_all,1);
xlim(ax, vc(1) + [-0.4*daxiao 0.4*daxiao]);
ylim(ax, vc(2) + [-0.4*daxiao 0.4*daxiao]);
zlim(ax, vc(3) + [-0.4*daxiao 0.4*daxiao]);
% xlim(ax, vc(1) + [-0.4*0.2 0.4*0.2]);
% ylim(ax, vc(2) + [-0.4*0.2 0.4*0.2]);
% zlim(ax, vc(3) + [-0.4*0.2 0.4*0.2]);
drawnow;

%% ==================== 动画主循环（同时对胸腔与心脏做形变，并进行雷达仿真） ====================
for frame = 1:total_frames
    t = (frame - 1) / fps;
    % 胸腔 M_t (呼吸四段式) — 与原胸腔代码完全一致
    t_mod_chest = mod(t, T_h_chest);
    if t_mod_chest < T_in
        M_t_chest = sin( (pi/2) * (t_mod_chest / T_in) )^2;
    elseif t_mod_chest < T_in + T_plateau
        M_t_chest = 1;
    elseif t_mod_chest < T_in + T_plateau + T_out
        x = (t_mod_chest - T_in - T_plateau) / T_out;
        M_t_chest = sin( (pi/2) + (pi/2) * x )^2;
    else
        M_t_chest = 0;
    end

    % 心脏 M_t (高斯尖峰) — 与原心脏代码一致
    t_mod_heart = mod(t, T_h_heart);
    M_t_heart = exp(-((t_mod_heart - tau_h_heart) / sigma_h_heart)^2);

    % 当前帧顶点位置（GPU）
    % 注意：delta_max_all_gpu 已经将 chest 和 heart 的 delta_max 拼接好，
    % 这里用一个系数向量把 chest 的系数设为 M_t_chest, heart 的系数设为 M_t_heart
    coeffs_gpu = [gpuArray(M_t_chest * ones(numVerts_chest,1)); gpuArray(M_t_heart * ones(numVerts_heart,1))];
    V_frame_gpu = Vertices_all_gpu + (coeffs_gpu .* delta_max_all_gpu);

    % 更新动画（取回CPU用于绘图）
    V_frame = gather(V_frame_gpu);
    set(p, 'Vertices', V_frame);
    title(ax, sprintf('Frame %d / %d (M_chest=%.3f, M_heart=%.3f)', frame, total_frames, M_t_chest, M_t_heart));
    drawnow;

    % ---------------------- 计算顶点速度 v_obj (基于帧间位移，GPU) ----------------------
    
    dt_frame = 1 / fps;  % 帧间时间间隔 [s]
    v_obj_gpu = (V_frame_gpu - V_prev_gpu) / dt_frame;  % Nx3 速度向量 [m/s]
    V_prev_gpu = V_frame_gpu;  % 更新上一帧

    % ---------------------- 雷达几何与反射计算（GPU） ----------------------

        d_i_gpu = V_frame_gpu - p_radar_gpu;                 % Nx3
        d_i_norm_gpu = sqrt(sum(d_i_gpu.^2, 2));            % r_i [m]
        d_hat_gpu = d_i_gpu ./ (d_i_norm_gpu + epsilon);     % 单位向量
    
        r_i_gpu = d_i_norm_gpu;                              % 距离 [m]
    
        v_i_gpu = sum(d_i_gpu .* v_obj_gpu, 2) ./ (d_i_norm_gpu + epsilon);  % 径向速度 [m/s]
    
        theta_i_gpu = atan2(d_i_gpu(:,2), d_i_gpu(:,1)) - theta_0;  % 方位角 [rad]
        phi_i_gpu   = asin(d_i_gpu(:,3) ./ (r_i_gpu + epsilon)) - phi_0;  % 俯仰角 [rad]
    
        % Lambertian 反射系数：接收增益（GPU）
        s_i_gpu = max(0, -sum(d_hat_gpu .* vertex_normals_all_gpu, 2));  % 0~1
    
    
        % ==================== 新增：胸腔与心脏反射能量加权 ====================
        w_chest = 1;   % 胸腔反射权重（可调）
        w_heart = 0.01;   % 心脏反射权重（可调）
    
        s_weights = [w_chest * ones(numVerts_chest,1,'gpuArray'); ...
                     w_heart * ones(numVerts_heart,1,'gpuArray')];
    
        s_i_gpu = s_i_gpu .* s_weights;
        % ================================================================
    
    
        % ---------------------- FMCW信号生成（多目标，解析信号版本，GPU加速） ----------------------
        mix_frame_gpu = zeros(1, length(t_frame_fmcw_gpu), 'gpuArray');  % 复数基带信号
    
        numSamples = length(t_frame_fmcw_gpu);
    
        r0_point_gpu = r_i_gpu;  % numVerts x 1
        v_point_gpu = v_i_gpu;   % numVerts x 1
    
        % 广播计算 r_t (numVerts x numSamples)
        r_t_gpu = r0_point_gpu + v_point_gpu .* t_frame_fmcw_gpu;  % 广播
        td_gpu = 2 * r_t_gpu / c;                                  % numVerts x numSamples
    
        % 发射相位 (1 x numSamples)
        phi_tx_gpu = 2 * pi * (fc * t_frame_fmcw_gpu + 0.5 * slope * t_frame_fmcw_gpu.^2);
    
        % 接收相位 (numVerts x numSamples)
        phi_rx_gpu = 2 * pi * (fc * (t_frame_fmcw_gpu - td_gpu) + 0.5 * slope * (t_frame_fmcw_gpu - td_gpu).^2);
    
        % 解析信号
        Tx_gpu = exp(1j * phi_tx_gpu);  % 1 x numSamples
        Rx_gpu = exp(1j * phi_rx_gpu) .* s_i_gpu;  % numVerts x numSamples (广播 s_i_gpu)
    
        % 混频 (numVerts x numSamples)
        mix_point_gpu = Tx_gpu .* conj(Rx_gpu);  % 广播 Tx_gpu
    
        % 累加所有点 (1 x numSamples)
        mix_frame_gpu = sum(mix_point_gpu, 1);
    
        % 转成 Nsample × Nchirp 并写入本帧 (取回CPU)
        mix_frame = gather(mix_frame_gpu);
        Mix(:,:,frame) = reshape(mix_frame, Nsample, Nchirp);

end

fprintf('✅ 动画与 FMCW 雷达仿真完成: %d 帧, %d FPS\n', total_frames, fps);
fprintf('Mix 信号矩阵大小: (Nsample=%d × Nchirp=%d × Nframe=%d)\n', Nsample, Nchirp, Nframe);

%% ==================== 后处理与可视化（保留你原来的 FFT / RD 展示逻辑） ====================
% 第一个 chirp 的距离 FFT（与胸腔代码一致）
Mix_gpu = gpuArray(Mix(:,1,:));
Mix_gpu = squeeze(Mix_gpu);
sig_fft1_gpu = fft(Mix_gpu, [], 1);   % 对第1维做FFT
sig_fft = gather(abs(sig_fft1_gpu));  % 取回CPU

figure;
mesh((0:total_frames-1)/1,(0:Nsample-1)*(maxR/Nsample),squeeze((sig_fft)));
ylabel('距离（频率bin）');
xlabel('时间');
zlabel('幅度');
title('距离维FFT结果 (第一个chirp)');

% 第一个帧的距离-多普勒图（与胸腔代码一致）
Mix_gpu = gpuArray(Mix(:,:,1));  % 第1帧
Mix_gpu = squeeze(Mix_gpu);
sig_fft1_gpu = fft(Mix_gpu, [], 1);   % 距离向FFT
sig_fft1_gpu = fftshift(fft(sig_fft1_gpu, [], 2),2); % 多普勒
sig_fft = gather(abs(sig_fft1_gpu));  % 取回CPU
figure;
mesh(squeeze(sig_fft));
ylabel('距离（频率bin）');
xlabel('速度');
zlabel('幅度');
title('距离-多普勒 (第一帧)');

%{
Mix_gpu = gpuArray(Mix(:,1,:));
Mix_gpu = squeeze(Mix_gpu);
% % --- 加窗（Hamming窗） ---
% win = gpuArray(hamming(Nsample));     % 列向量窗
% Mix_gpu = Mix_gpu .* win;             % 窗函数作用在快时间维

sig_fft1_gpu = fft(Mix_gpu, [], 1);   % 对第1维做FFT
sig_fft = gather((sig_fft1_gpu));  % 取回CPU
% sig_fft = squeeze(sig_fft) - mean(squeeze(sig_fft),2);
RangeFFTcor = squeeze(sig_fft);
[~, MaxIndex2] = max(mean(RangeFFTcor,2));
It2 = zeros(1,total_frames);
Qt2 = zeros(1,total_frames);

for ii = 1: total_frames

    It2(ii) = real(RangeFFTcor(MaxIndex2,ii));
    Qt2(ii) = imag(RangeFFTcor(MaxIndex2,ii));
end


N = total_frames;
for n=1:N
    phi2(n)=0;
    for k=2:n
        phi2(n)=phi2(n)+(It2(k)*(Qt2(k)-Qt2(k-1))-Qt2(k)*(It2(k)-It2(k-1)))/(It2(k).^2+Qt2(k).^2);
    end
end

figure;
plot(phi2);


angleDenoised=diff(phi2);


%plot phase data and frequency spectrum
figure;
subplot(121)
plot((0:length(angleDenoised)-1)*(total_seconds/Nframe),angleDenoised,'LineWidth',2)
xlabel("时间 (s)");
ylabel("相位 (rads)");
title("差分相位序列时序图");

figure;
plot((-length(angleDenoised)/2:length(angleDenoised)/2-1)*fps/length(angleDenoised),abs(fftshift(fft(angleDenoised)))/max(abs(fftshift(fft(angleDenoised)))),'LineWidth',2)
axis([0 10 0 1.1]);
xlabel("频率 (Hz)");
ylabel("幅度");
title("差分相位序列频率图");
%}

