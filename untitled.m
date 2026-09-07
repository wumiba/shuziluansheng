clear; close all; clc;

%% ---------------------- STL 文件名（按需修改） ----------------------
chest_stl = 'd.stl';   % 胸腔 STL（若你原来用 d.stl 可改回）
heart_stl = 'b_三倍大小.stl';   % 心脏 STL（若你原来用 b.stl 可改回）

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
d_h_chest      = 0.07;       % 最大形变位移量 [m]
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
total_seconds = 4;
total_frames = round(total_seconds * fps);

% chest stl 路径（使用 chest_stl 变量）

use_vertex_normals_chest = true; % chest 原脚本 true

%% ==================== === 心脏（来源：你给出的“心脏”代码） ====================
% 前缀 heart_
daxiao_heart = 3;
T_h_heart    = 0.8;        % 心跳周期（秒）
fps_heart    = 25;         % 帧率（与胸腔保持一致）
tau_h_heart  = T_h_heart / 3;
sigma_h_heart = T_h_heart / 5;
d_h_heart     = 0.003*daxiao_heart;     % 最大形变位移量 [m]
direction_heart = [-0.5, -1.0, 0.0];

% 心脏影响范围（心脏脚本里的）
x_min_heart = -0.016; x_max_heart = 0.03;
y_min_heart = -0.02;  y_max_heart = 0.01;
z_min_heart = -0.022; z_max_heart = 0.0075;

x_bian = 0.10;
z_bian = 0.14;

x_min_heart = x_min_heart*daxiao_heart+x_bian; x_max_heart = x_max_heart*daxiao_heart+x_bian;
y_min_heart = y_min_heart*daxiao_heart;  y_max_heart = y_max_heart*daxiao_heart;
z_min_heart = z_min_heart*daxiao_heart+z_bian; z_max_heart = z_max_heart*daxiao_heart+z_bian;
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
Mix = zeros(N_rx,Nsample, Nchirp, Nframe);

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
Vertices_heart = Vertices_heart+[x_bian 0 z_bian];

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
%% ==================== 绘图初始化 ====================
fig = figure('Name','Chest + Heart Animation with Radar Echo','Color',[1 1 1]);
ax = axes('Parent', fig);
hold(ax, 'on');

% --- 胸腔 patch（蓝色半透明） ---
p_chest = patch('Vertices', Vertices_chest, 'Faces', Faces_chest, ...
    'FaceColor',[0.3 0.5 1.0], ...       % 蓝色
    'FaceAlpha',0.5, ...                % 半透明
    'EdgeColor','none', ...
    'FaceLighting','gouraud', ...
    'AmbientStrength',0.3, ...
    'SpecularStrength',0.2, ...
    'Parent', ax);

% --- 心脏 patch（红色不透明） ---
p_heart = patch('Vertices', Vertices_heart, 'Faces', Faces_heart, ...
    'FaceColor',[0.8 0.2 0.2], ...
    'FaceAlpha',1.0, ...
    'EdgeColor','none', ...
    'FaceLighting','gouraud', ...
    'AmbientStrength',0.3, ...
    'SpecularStrength',0.2, ...
    'Parent', ax);

axis equal off;                % 保持比例但关闭坐标轴显示
axis vis3d;                    % 防止视图自动缩放
material shiny;                % 保留原材质
camlight headlight;            % 保留光照方向
view(ax, [1 0 0]);             % 从 +X 方向看向 -X
hold(ax, 'off');







%% ==================== 导出 GIF（固定从 +X 看向 -X） ====================
gif_name = 'ppt_4s.gif';
fps_gif = 25;  % 与动画帧率保持一致
delay_time = 1 / fps_gif;

% 设置相机视角：从 +X 看向 -X
view(ax, [1 0 0]);
camlight(ax, 'headlight'); % 保留你原来的光照方向
drawnow;

% 重新播放并录制成 GIF（保持原显示风格）
for frame = 1:total_frames
    t = (frame - 1) / fps;

    % 重复使用上面的形变计算逻辑
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

    t_mod_heart = mod(t, T_h_heart);
    M_t_heart = exp(-((t_mod_heart - tau_h_heart) / sigma_h_heart)^2);

    coeffs_gpu = [gpuArray(M_t_chest * ones(numVerts_chest,1)); gpuArray(M_t_heart * ones(numVerts_heart,1))];
    V_frame_gpu = Vertices_all_gpu + (coeffs_gpu .* delta_max_all_gpu);
    V_frame = gather(V_frame_gpu);

    V_chest = V_frame(1:numVerts_chest, :);
    V_heart = V_frame(numVerts_chest+1:end, :);
    set(p_chest, 'Vertices', V_chest);
    set(p_heart, 'Vertices', V_heart);
    xlim([-0.2, 0.2]);
    ylim([-0.4, 0.4]);
    zlim([-0.6, 0.7]);

    drawnow;

    frame_img = getframe(fig);
    [A, map] = rgb2ind(frame2im(frame_img), 256);

    if frame == 1
        imwrite(A, map, gif_name, 'gif', 'LoopCount', inf, 'DelayTime', delay_time);
    else
        imwrite(A, map, gif_name, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end
end

disp(['✅ GIF 已保存为 ', gif_name]);
