
clear; close all; clc;

% %% ---------------------- 心跳用户参数 ----------------------
% T_h      = 2;        % 心跳周期（秒）
% fps      = 25;         % 帧率（动画帧率，与FMCW帧对应）
% tau_h    = T_h / 3;    % 心跳强度最大时刻
% sigma_h  = T_h / 5;   % 控制心跳尖锐程度
% d_h      = 0.03;      % 最大形变位移量 [m]
% direction = [-0.5, -1.0, 0.0];  % 固定形变方向向量

%% ---------------------- 呼吸用户参数 ----------------------
% —— 呼吸分段参数（按图）
T_in       = 3;      % 吸气上升阶段时长 [s]
T_plateau  = 0.8;      % 平台期（吸气或呼气的维持）[s]
T_out      = 2;      % 呼气下降阶段时长 [s]
T_b        = T_in + T_plateau + T_out + T_plateau;   % 完整呼吸周期
T_h      = T_b;        % 用呼吸周期覆盖原心跳周期（其余代码仍用 T_h）
fps      = 25;         % 帧率
% tau_h、sigma_h 行可以保留但后续不会再用
tau_h    = T_h / 3;
sigma_h  = T_h / 5;
d_h      = 0.03;       % 最大形变位移量 [m]
direction = [-0.5, -1.0, 0.0];


% 心脏影响范围
x_min = -0.13; x_max = 0.13;
y_min = -0.15; y_max = 0.08;
z_min = -0.2; z_max = 0.2;

% 动画长度
total_seconds = 15;
total_frames = round(total_seconds * fps);  % 总帧数（Nframe = total_frames）

% STL 文件路径
stl_file = 'd.stl';

% 是否使用顶点法向量
use_vertex_normals = true; % true以法向量形式扩张，false沿固定方向扩张

%% ---------------------- 毫米波FMCW雷达参数设置 ----------------------
fc = 77e9;             % 雷达工作频率 [Hz]
c = 3e8;               % 光速 [m/s]

maxR = 200;            % 最大探测距离 [m]
rangeRes = 1;          % 距离分辨率 [m]
maxV = 20;             % 最大速度 [m/s]

B = c / (2 * rangeRes); % 带宽 [Hz]
Nchirp = 8;              % chirp数（每个FMCW帧内的chirp数）
Nsample = 128;              % 每个chirp的采样点数
Nframe = total_frames;     % 帧数（与动画帧对应）

Tchirp = 2.5 * 2 * maxR / c;  % chirp持续时间 [s]
idle_time = 5e-6;    % idle时间 [s]
vres = (c / fc) / (2 * Nchirp * (Tchirp + idle_time)); % 速度分辨率 [m/s]

slope = B / Tchirp;    % 调频斜率 [Hz/s]

% 雷达位置和指向（假设雷达位于心脏前方）
p_radar = [0, 5, 0]; % 雷达位置 [m]（x,y,z）
theta_0 = 0;           % 雷达方位角参考 [rad]
phi_0 = 0;             % 雷达俯仰角参考 [rad]

epsilon = 1e-12;       % 数值稳定性

% 数据存储：三维FMCW混频信号 (Nsample × Nchirp × Nframe)
Mix = zeros(Nsample, Nchirp, Nframe);  % 复数基带信号（最终从GPU取回）

%% ---------------------- 初始化GPU ----------------------
% 检查并选择GPU设备
if gpuDeviceCount > 0
    gpu = gpuDevice(1);  % 使用第一个可用GPU
    disp('✅ 使用GPU加速: ' + string(gpu.Name));
else
    error('❌ 未检测到兼容GPU。请检查硬件和MATLAB Parallel Computing Toolbox。');
end

%% ---------------------- 读取 STL ----------------------
[Vertices, Faces] = read_stl_generic(stl_file);
Vertices = double(Vertices);
Faces = double(Faces);



%% ---------------------- 计算顶点法向量 ----------------------
numFaces = size(Faces,1);
numVerts = size(Vertices,1);
n_sum = zeros(numVerts,3);
for f = 1:numFaces
    i1 = Faces(f,1); i2 = Faces(f,2); i3 = Faces(f,3);
    v1 = Vertices(i1,:); v2 = Vertices(i2,:); v3 = Vertices(i3,:);
    nf = cross(v2 - v1, v3 - v1);
    nf = nf / (norm(nf) + epsilon);
    n_sum(i1,:) = n_sum(i1,:) + nf;
    n_sum(i2,:) = n_sum(i2,:) + nf;
    n_sum(i3,:) = n_sum(i3,:) + nf;
end
vertex_normals = n_sum ./ (sqrt(sum(n_sum.^2,2)) + epsilon); % Nx3

% 将顶点法向量移到GPU（后续使用）
vertex_normals_gpu = gpuArray(vertex_normals);

%% ---------------------- 计算形变幅度 A_i ----------------------
e_const = exp(1);
A = zeros(numVerts,1);
for i=1:numVerts
    x = Vertices(i,1); y = Vertices(i,2); z = Vertices(i,3);
    if ~(x >= x_min && x <= x_max && y >= y_min && y <= y_max && z >= z_min && z <= z_max)
        continue;
    end
    fz = - ( (z - z_min) * (z - z_max) ) / ((z_max - z_min)^2);
    fz = max(fz, 0);
    fy = - ( (y - y_min) * (y - y_max) ) / ((y_max - y_min)^2);
    fy = max(fy, 0);
    sqrt_fz = fz;
    sqrt_fy = sqrt(fy);
    ratio = (x - x_min) / (x_max - x_min);
    ratio = min(max(ratio,0),1);
    log_term = log( e_const + (1 - e_const) * ratio );
    A(i) = d_h * sqrt_fz * 1 * log_term;
end

% 将A移到GPU
A_gpu = gpuArray(A);

%% ---------------------- 确定形变方向 ----------------------
if use_vertex_normals
    delta_max_gpu = A_gpu .* vertex_normals_gpu;  % GPU上计算
else
    direction = direction / norm(direction); % 归一化
    direction_gpu = gpuArray(direction);
    delta_max_gpu = A_gpu .* direction_gpu;
end

%% ---------------------- 绘图初始化 ----------------------
fig = figure('Name','Heart Beat Animation with Radar Echo','Color',[1 1 1]);
ax = axes('Parent', fig);
p = patch('Vertices', Vertices, 'Faces', Faces, ...
    'FaceColor',[0.8 0.3 0.3], 'EdgeColor','none', ...
    'FaceLighting','gouraud', 'AmbientStrength',0.3, ...
    'SpecularStrength',0.2, 'Parent', ax);
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
material shiny; camlight headlight; view(3); grid on;
vc = mean(Vertices,1);
xlim(ax, vc(1) + [-0.4 0.4]);
ylim(ax, vc(2) + [-0.4 0.4]);
zlim(ax, vc(3) + [-0.4 0.4]);
% axis auto;
drawnow;

%% ---------------------- 动画与FMCW雷达仿真（GPU加速） ----------------------
% 初始化上一帧顶点位置（GPU）
Vertices_gpu = gpuArray(Vertices);
V_prev_gpu = Vertices_gpu;
p_radar_gpu = gpuArray(p_radar);

% FMCW时间轴（GPU）
t_frame_fmcw_gpu = gpuArray(linspace(0, Nchirp * Tchirp, Nsample * Nchirp));  % [s]

for frame = 1:total_frames
    t = (frame - 1) / fps;  % 时间 [s]
    t_mod = mod(t, T_h);
%     M_t = exp(-((t_mod - tau_h) / sigma_h)^2);


    % —— 呼吸调制 M_t（按图的四段）
if t_mod < T_in
    % 上升段：sin^2( (pi/2) * (t_mod/T_in) )
    M_t = sin( (pi/2) * (t_mod / T_in) )^2;
elseif t_mod < T_in + T_plateau
    % 吸气平台期
    M_t = 1;
elseif t_mod < T_in + T_plateau + T_out
    % 下降段：sin^2( pi/2 + (pi/2) * ((t_mod - T_in - T_plateau)/T_out) )
    x = (t_mod - T_in - T_plateau) / T_out;
    M_t = sin( (pi/2) + (pi/2) * x )^2;
else
    % 呼气平台期
    M_t = 0;
end

    
    % 当前帧顶点位置（GPU）
    V_frame_gpu = Vertices_gpu + (M_t .* delta_max_gpu);


    %test
    %{
figure;
    for frame = 1:total_frames
    t = (frame - 1) / fps;  % 时间 [s]
    t_mod = mod(t, T_h);
%     M_t = exp(-((t_mod - tau_h) / sigma_h)^2);


    % —— 呼吸调制 M_t（按图的四段）
if t_mod < T_in
    % 上升段：sin^2( (pi/2) * (t_mod/T_in) )
    M_t = sin( (pi/2) * (t_mod / T_in) )^2;
elseif t_mod < T_in + T_plateau
    % 吸气平台期
    M_t = 1;
elseif t_mod < T_in + T_plateau + T_out
    % 下降段：sin^2( pi/2 + (pi/2) * ((t_mod - T_in - T_plateau)/T_out) )
    x = (t_mod - T_in - T_plateau) / T_out;
    M_t = sin( (pi/2) + (pi/2) * x )^2;
else
    % 呼气平台期
    M_t = 0;
end
        aaa(frame) = M_t;
    end
    plot(aaa);
    
    
    %}

    
    % 更新动画（取回CPU用于绘图）
    V_frame = gather(V_frame_gpu);  % 从GPU取回
    set(p, 'Vertices', V_frame);
    title(ax, sprintf('Heart Apex Beating — frame %d / %d (M_t=%.3f)', ...
        frame, total_frames, M_t));
    drawnow;
%     pause(1/fps);
    
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
    s_i_gpu = max(0, -sum(d_hat_gpu .* vertex_normals_gpu, 2));  % 0~1
    
    % ---------------------- FMCW信号生成（多目标，解析信号版本，GPU加速） ----------------------
    mix_frame_gpu = zeros(1, length(t_frame_fmcw_gpu), 'gpuArray');  % 复数基带信号
    
    % 向量化所有顶点（避免循环，使用广播）
    % 扩展维度：r0_point_gpu (numVerts x 1), v_point_gpu (numVerts x 1), s_i_gpu (numVerts x 1)
    % t_frame_fmcw_gpu (1 x numSamples)
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

fprintf('✅ 动画完成: %d 帧, %d FPS\n', total_frames, fps);
fprintf('✅ FMCW雷达仿真完成: Mix信号矩阵 (Nsample=%d × Nchirp=%d × Nframe=%d)\n', Nsample, Nchirp, Nframe);

%% ---------------------- 第一个chirp的距离FFT ----------------------
Mix_gpu = gpuArray(Mix(:,1,:));  % 第个chirp
% % --- 加窗（Hamming窗） ---
% win = gpuArray(hamming(Nsample));     % 列向量窗
% Mix_gpu = Mix_gpu .* win;             % 窗函数作用在快时间维

sig_fft1_gpu = fft(Mix_gpu, [], 1);   % 对第1维做FFT
sig_fft = gather(abs(sig_fft1_gpu));  % 取回CPU

sig_fft = squeeze(sig_fft) - mean(squeeze(sig_fft),2);

figure;
mesh((0:total_frames-1)/1,(0:Nsample-1)*(maxR/Nsample),squeeze((sig_fft)));
ylabel('距离（频率bin）');
xlabel('时间');
zlabel('幅度');
title('距离维FFT结果 (第一个chirp)');



%% ---------------------- 第一个chirp帧的距离多普勒---------
Mix_gpu = gpuArray(Mix(:,:,1));  % 第个chirp
sig_fft1_gpu = fft(Mix_gpu, [], 1);   % 对第1维做FFT
sig_fft1_gpu = fftshift(fft(sig_fft1_gpu, [], 2),2);
sig_fft = gather(abs(sig_fft1_gpu));  % 取回CPU
figure;
mesh(squeeze(sig_fft));
ylabel('距离（频率bin）');
xlabel('速度');
zlabel('幅度');
title('距离维FFT结果 (第一帧)');



% Mix_gpu = gpuArray(Mix(:,:,1));              % [Nsample x Nchirp]
% 
% % 窗函数
% w_r = gpuArray(hann(size(Mix_gpu,1)));       % 距离维窗
% w_d = gpuArray(hann(size(Mix_gpu,2))).';     % 多普勒维窗
% X = Mix_gpu .* w_r .* w_d;                   % 双向加窗（按需）
% 
% % 距离维 FFT
% RD = fft(X, [], 1);
% 
% % 多普勒维 FFT
% RD = fftshift(fft(RD, [], 2), 2);
% 
% RD_abs = gather((abs(RD)+eps));
% figure; mesh(RD_abs); axis xy
% xlabel('速度 bin'); ylabel('距离 bin');
% title('Range–Doppler Map'); colorbar



%% ---------------------- STL 读取函数 ----------------------
function [V, F] = read_stl_generic(filename)
    fid = fopen(filename,'r');
    if fid < 0, error('Cannot open file: %s', filename); end
    header = fread(fid,80,'uint8=>char')';
    frewind(fid);
    is_ascii = startsWith(strtrim(lower(header)), 'solid');
    if is_ascii
        fclose(fid);
        str = fileread(filename);
        vertexTokens = regexp(str, 'vertex\s+([-+eE0-9\.-]+)\s+([-+eE0-9\.-]+)\s+([-+eE0-9\.-]+)', 'tokens');
        verts = cell2mat(cellfun(@(c) str2double(c), vertexTokens, 'UniformOutput', false));
        numVerts = size(verts,1);
        if mod(numVerts,3) ~= 0, error('顶点数不是3的倍数'); end
        rawV = verts;
        rawF = reshape(1:numVerts, 3, [])';
    else
        fseek(fid,80,'bof');
        numFaces = fread(fid,1,'uint32');
        rawV = zeros(numFaces*3,3);
        rawF = reshape(1:(numFaces*3),3,[])';
        for k=1:numFaces
            fread(fid,3,'float32'); % normal
            v1 = fread(fid,3,'float32')';
            v2 = fread(fid,3,'float32')';
            v3 = fread(fid,3,'float32')';
            rawV((k-1)*3+1,:) = v1;
            rawV((k-1)*3+2,:) = v2;
            rawV((k-1)*3+3,:) = v3;
            fread(fid,1,'uint16'); % attr
        end
        fclose(fid);
    end
    tol = 1e-9;
    [uniqueV,~,idxMap] = uniquetol(rawV, tol, 'ByRows', true);
    V = uniqueV; F = reshape(idxMap, 3, [])';
end