% fused_heartbeat_radar_simulation_gpu_stl.m
% 心脏 STL 形变动画 + exp 解析信号 FMCW 雷达仿真 + GPU 加速（批处理）
% 输出：动画 + MixVolume (Nr x Nd x Nf)
% 图像：1) 第一帧第一个chirp距离FFT；2) 第一帧距离FFT谱；
%       3) 第一帧距离-多普勒谱；4) 第一个chirp 距离×时间（帧）

clear; close all; clc;

%% ===================== GPU 设置 =====================
gpu_enable = true;           % 有GPU就开，没有会自动回退
gpu_batch_vertices = 4000;   % 每批顶点数（按显存调）
template = complex(single(0), single(0));  % 复数模板：解决 gpuArray.zeros 维度报错

if gpu_enable
    try
        g = gpuDevice;
        fprintf('✅ 使用GPU设备: %s (CC %s, %.1f GB)\n', g.Name, g.ComputeCapability, g.TotalMemory/1024^3);
    catch
        warning('⚠️ 未检测到可用GPU，自动回退CPU。');
        gpu_enable = false;
    end
end

rng(1); % 固定随机种子

%% ===================== 心脏动画参数（来自代码1） =====================
T_h      = 0.8;        % 心跳周期（秒）
fps      = 25;         % 帧率
tau_h    = T_h / 3;    % 心跳强度最大时刻
sigma_h  = T_h / 10;   % 控制心跳尖锐程度

d_h      = 0.003;      % 最大形变位移量
direction = [-0.5, -1.0, 0.0];  % 固定形变方向向量

% 心脏影响范围
x_min = -0.016; x_max = 0.03;
y_min = -0.02;  y_max = 0.01;
z_min = -0.022; z_max = 0.0075;

% 动画长度
total_seconds = 15;
total_frames  = round(total_seconds * fps);

% STL 文件路径
stl_file = 'b.stl';

% 是否使用顶点法向量
use_vertex_normals = true;

%% ===================== 读取 STL（沿用代码1） =====================
[Vertices, Faces] = read_stl_generic(stl_file);
Vertices = double(Vertices);
Faces    = double(Faces);

%% ===================== 顶点法向量（沿用代码1） =====================
epsilon  = 1e-12;
numFaces = size(Faces,1);
numVerts = size(Vertices,1);
n_sum    = zeros(numVerts,3);

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

%% ===================== 形变幅度 A_i（沿用代码1） =====================
e_const = exp(1);
A = zeros(numVerts,1);
for i=1:numVerts
    x = Vertices(i,1); y = Vertices(i,2); z = Vertices(i,3);
    if ~(x >= x_min && x <= x_max && y >= y_min && y <= y_max && z >= z_min && z <= z_max), continue; end
    fz = - ( (z - z_min) * (z - z_max) ) / ((z_max - z_min)^2);
    fz = max(fz, 0);
    fy = - ( (y - y_min) * (y - y_max) ) / ((y_max - y_min)^2);
    fy = max(fy, 0);
    sqrt_fz = sqrt(fz);
    sqrt_fy = sqrt(fy);
    ratio = (x - x_min) / (x_max - x_min);
    ratio = min(max(ratio,0),1);
    log_term = log( e_const + (1 - e_const) * ratio );
    A(i) = d_h * sqrt_fz * sqrt_fy * log_term;
end

%% ===================== 形变方向（沿用代码1） =====================
if use_vertex_normals
    delta_max = A .* vertex_normals;
else
    direction = direction / norm(direction);
    delta_max = A .* direction;
end

%% ===================== 雷达系统参数（exp 解析信号） =====================
fc = 77e9;         % 载频
c  = 3e8;
rangeRes = 1;      % m
B   = c/(2*rangeRes);
Nd  = 128;         % chirp 数
Nr  = 256;         % 每chirp采样点
Nf  = total_frames;% 帧数（与动画同步）
Tchirp = 5 * 2 * 200 / c;  % 基于 maxR=200 m
slope  = B / Tchirp;
Fs     = Nr / Tchirp;
endle_time = 2.5e-6;
vres   = (c/fc)/(2*Nd*(Tchirp+endle_time));

% 单帧时间轴（Nr*Nd 个采样）
Nt = single(Nr * Nd);
t_frame_cpu = linspace(0, Nd*Tchirp, Nt);
if gpu_enable, t_frame = gpuArray(single(t_frame_cpu)); else, t_frame = single(t_frame_cpu); end

% 发射相位（所有帧相同）
phi_tx = 2*pi*( single(fc) .* t_frame + 0.5*single(slope) .* (t_frame.^2) );
Tx     = exp(1j * phi_tx);  % 1 x Nt (GPU/CPU)

% 雷达位置与指向（在心脏模型前方，指向质心）
vc = mean(Vertices,1);
p_radar = vc + [0.30, 0.0, 0.0];   % +30 cm 沿 x
vec_to_model = (vc - p_radar);
theta_0 = atan2(vec_to_model(2), vec_to_model(1));
phi_0   = asin(vec_to_model(3) / (norm(vec_to_model)+epsilon));

% 接收天线增益（线性）
G_rx = single(4);  % ~6 dBi
% 噪声
noise_sigma = single(1e-3);

%% ===================== 结果容器 =====================
% 替换后
%% ===================== 结果容器 =====================
if gpu_enable
    % 修复：先在 CPU 上创建 zeros 矩阵，再移动到 GPU
    MixVolume = gpuArray(zeros(Nr, Nd, Nf, 'like', template));
else
    MixVolume = zeros(Nr, Nd, Nf, 'like', template);
end

%% ===================== 动画初始化（可关） =====================
fig = figure('Name','Heart Beat + Radar Echo (GPU)','Color',[1 1 1]);
ax  = axes('Parent', fig);
p   = patch('Vertices', Vertices, 'Faces', Faces, ...
            'FaceColor',[0.8 0.3 0.3], 'EdgeColor','none', ...
            'FaceLighting','gouraud', 'AmbientStrength',0.3, ...
            'SpecularStrength',0.2, 'Parent', ax);
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
material shiny; camlight headlight; view(3); grid on;
xlim(ax, vc(1) + [-0.05 0.05]);
ylim(ax, vc(2) + [-0.05 0.05]);
zlim(ax, vc(3) + [-0.05 0.05]);
drawnow;

%% ===================== 逐帧：形变 + 回波合成 (GPU批处理) =====================
r_prev = zeros(numVerts,1,'single');  % 上一帧距离

for frame = 0:(total_frames-1)
    % 心跳调制
    t      = frame / fps;
    t_mod  = mod(t, T_h);
    M_t    = exp(-((t_mod - tau_h) / sigma_h)^2);
    % 顶点形变
    V_frame = Vertices + (M_t .* delta_max);
    % 动画
    set(p, 'Vertices', V_frame);
    title(ax, sprintf('Heart Beat — frame %d/%d (M_t=%.3f)', frame+1, total_frames, M_t));
    drawnow limitrate;

    % 几何
    d_i   = V_frame - p_radar;           % Nx3
    r_i   = sqrt(sum(d_i.^2,2));         % Nx1 (double)
    d_hat = d_i ./ (r_i + epsilon);      % Nx3

    % Lambertian + 接收增益 + 1/r^2
    s_i = max(0, -sum(d_hat .* vertex_normals, 2));   % Nx1
    alpha_k = (s_i .* double(G_rx)) ./ ((r_i + epsilon).^2); % double

    % 转成 single，准备上 GPU
    r_i_s   = single(r_i);
    alpha_s = single(alpha_k);

    % GPU/CPU 回波叠加
    if gpu_enable
        Rx_total = gpuArray(zeros(1, Nt, 'like', template));  % 1 x Nt (复)
        % 按顶点批处理（避免 K×Nt 巨矩阵）
        for k0 = 1:gpu_batch_vertices:numVerts
            k1 = min(k0 + gpu_batch_vertices - 1, numVerts);
            td_k = 2 * r_i_s(k0:k1) / single(c);  % Kx1
            a_k  = alpha_s(k0:k1);                % Kx1

            % K x Nt：利用隐式扩展
            tau   = t_frame - td_k;
            valid = tau >= 0;
            % 相位
            phi_rx = 2*pi*( single(fc).*tau + 0.5*single(slope).*(tau.^2) );
            contrib = a_k .* exp(1j*phi_rx);
            contrib(~valid) = 0;

            % 沿顶点维求和 -> 1 x Nt
            Rx_batch = sum(contrib, 1);
            Rx_total = Rx_total + Rx_batch;
        end

        % 噪声
        Rx_total = Rx_total + (noise_sigma/sqrt(2))*(randn(size(Rx_total),'like',template) + 1j*randn(size(Rx_total),'like',template));
        % 混频：exp解析信号，共轭解调
        mix_frame = Tx .* conj(Rx_total); % 1 x Nt (GPU)
    else
        Rx_total = zeros(1, Nt, 'like', template);
        for k0 = 1:gpu_batch_vertices:numVerts
            k1 = min(k0 + gpu_batch_vertices - 1, numVerts);
            td_k = 2 * r_i_s(k0:k1) / single(c);
            a_k  = alpha_s(k0:k1);

            tau   = single(t_frame_cpu) - td_k;     % 使用 CPU 轴
            valid = tau >= 0;
            phi_rx = 2*pi*( single(fc).*tau + 0.5*single(slope).*(tau.^2) );
            contrib = a_k .* exp(1j*phi_rx);
            contrib(~valid) = 0;

            Rx_total = Rx_total + sum(contrib, 1);
        end
        Rx_total = Rx_total + (noise_sigma/sqrt(2))*(randn(size(Rx_total),'like',template) + 1j*randn(size(Rx_total),'like',template));
        Tx_cpu   = exp(1j * gather(phi_tx));   % 防止重复计算
        mix_frame = Tx_cpu .* conj(Rx_total);  % CPU
    end

    % 存储 Nr x Nd
    MixVolume(:,:,frame+1) = reshape(mix_frame, Nr, Nd);

    % 控制播放速度（想更快可注释）
    pause(1/fps);
end

if gpu_enable
    fprintf('✅ 仿真完成（GPU）: MixVolume = [%d, %d, %d]\n', size(MixVolume));
else
    fprintf('✅ 仿真完成（CPU）: MixVolume = [%d, %d, %d]\n', size(MixVolume));
end

%% ===================== 后处理：只分析一部分帧 =====================
frames_to_analyze = 1 : min(16, Nf);  % 可改
Nf_an = numel(frames_to_analyze);

% 1) 第一帧第一个 chirp 的 距离FFT 曲线
first_frame = 1; first_chirp_idx = 1;
sig_firstchirp = MixVolume(:, first_chirp_idx, first_frame);
range_fft_firstchirp = abs( fft(gather(sig_firstchirp)) );

figure;
plot(range_fft_firstchirp);
xlabel('距离（频率bin）'); ylabel('幅度');
title('第一个chirp的FFT结果 (第一帧)');

% 2) 距离维FFT结果谱矩阵（第一帧）
sig_fft1 = fft(gather(MixVolume(:,:,first_frame)), [], 1); % Nr x Nd
sig_fft  = abs(sig_fft1);
figure;
mesh(sig_fft);
ylabel('距离（频率bin)'); xlabel('chirp脉冲数'); zlabel('幅度');
title('距离维FFT结果 (第一帧)');

% 3) 速度维FFT（第一帧） -> 距离-多普勒谱
sig_fft2 = fft(sig_fft1, [], 2); % 对 chirp 做 FFT
RDM = abs(sig_fft2);
doppler_axis = linspace(0,128,Nd)*vres;
range_axis   = linspace(0,256,Nr)*rangeRes;
figure;
mesh(doppler_axis, range_axis, RDM);
xlabel('多普勒通道'); ylabel('距离通道'); zlabel('幅度');
title('速度维FFT 距离-多普勒谱 (第一帧)');

% 4) 新增：第一个chirp 的 距离 × 时间（帧）
range_time_map = zeros(Nr, Nf_an);
for idx = 1:Nf_an
    f = frames_to_analyze(idx);
    ch = abs( fft(gather(MixVolume(:, first_chirp_idx, f))) );
    range_time_map(:, idx) = ch;
end

figure;
imagesc(frames_to_analyze, range_axis, range_time_map);
set(gca,'YDir','normal');
xlabel('帧索引'); ylabel('距离 (m)');
title(sprintf('第一个chirp 的 距离 × 时间（帧） 图 (前 %d 帧)', Nf_an));
colorbar;

%% ===================== STL 读取函数（与代码1相同） =====================
function [V, F] = read_stl_generic(filename)
    fid = fopen(filename,'r');
    if fid < 0, error('Cannot open file: %s', filename); end
    header = fread(fid,80,'uint8=>char')';
    frewind(fid);
    is_ascii = startsWith(strtrim(lower(header)), 'solid');
    if is_ascii
        fclose(fid);
        str = fileread(filename);
        vertexTokens = regexp(str, 'vertex\s+([-+eE0-9\.\-]+)\s+([-+eE0-9\.\-]+)\s+([-+eE0-9\.\-]+)', 'tokens');
        verts = cell2mat(cellfun(@(c) str2double(c), vertexTokens, 'UniformOutput', false));
        numVerts = size(verts,1);
        if mod(numVerts,3) ~= 0, error('顶点数不是3的倍数'); end
        rawV = verts; rawF = reshape(1:numVerts, 3, [])';
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
