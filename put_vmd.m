% put_vmd.m
% =========================================================================
% 自包含脚本 = wanzheng.m 的数字孪生 FMCW 回波仿真 + main_test.m 的 VMD 解耦出图
%
% 背景：把"胸腔呼吸 + 心脏搏动"两个 STL 三维模型做几何形变、合成 77GHz FMCW
%       多通道基带回波（与 wanzheng.m 完全一致），再从回波里经 Range FFT ->
%       相位解调 -> VMD(变分模态分解) -> 自动挑选呼吸/心搏 IMF，参照
%       main_test.m 的图与注释风格出图并保存 PNG。
%
% 用法：在 MATLAB 直接运行本文件即可（需 GPU + Parallel Computing Toolbox，
%       无需先运行 wanzheng.m）。STL 用默认 d.stl / b.stl（放同目录）。
%
% 说明1：仅保留几何形变与回波合成、去掉逐帧动画刷新（750 帧逐帧 patch 太慢，
%        且动画帧报告已有 fig5）。
% 说明2：VMD 需足够频率分辨率区分呼吸(~0.42 Hz)与心搏(~1.25 Hz)。帧数取
%        1000(=40 s x 25 fps)还有一个原因：VMD.m 的返回 imf 会把模态截取到
%        alpha=1000 列，输入不足 1000 点会越界，取 1000 帧可让输出为整段、
%        时域长度与频谱对齐（频率分辨率 0.025 Hz）。想更快可把 total_seconds 调小，
%        但需同步把 alpha 改成 <= 输入长度，避免 imf 截断报错。
% 说明3：main_test 里实测数据才需要的 f_phaseDenoise 脉冲去噪此处省略，
%        仿真回波无随机噪声；其调用方式以注释保留在相位提取处。
%
% 输出：PNG（300 dpi）存到  数字孪生素材_开题报告/图/（默认加噪 SNR=50 dB，σ_n=RMS/10^(SNR/20)）
%   组合图（保留）：fig6 全IMF时域、fig7 全IMF频谱、fig8 呼吸/心搏 2x2
%   fig9 谱峰对照（总体系一归一化，保留心搏/呼吸实际幅度比）
%   独立图（新增）：vmd_imf1..7_时域 / _频谱（每模态时域、频谱各一张）
%                vmd_呼吸_时域 / vmd_呼吸_频谱 / vmd_心跳_时域 / vmd_心跳_频谱
%                fig10 雷达距离-时间图（RX1 第一chirp，加噪后）
%   并在命令行打印 测得 BR/HR 与真值(60/T_b, 60/T_h_heart) 对照。
% =========================================================================

clear; close all; clc;

%% 目录与输出（素材包图目录）
thisdir  = fileparts(mfilename('fullpath'));
outdir   = fullfile(thisdir,'数字孪生素材_开题报告','图');
if ~exist(outdir,'dir'), mkdir(outdir); end

%% 仿真时长 / 慢时间采样率（fps 即 VMD 处理用到的 fs）
fps          = 25;          % 帧率 [Hz]（= 慢时间采样率 fs）
fs           = fps;         % 信号处理统一用 fs=fps
total_seconds= 40;          % 用于 VMD 的仿真时长 [s]（25 fps -> 1000 帧，见说明2）
total_frames = round(total_seconds*fps);
Nframe       = total_frames;

%% 数据来源开关：已用本脚本(含体表心搏脉冲)跑过并保存了无噪 Mix，
% 设 true 直接加载再按 SNR 加噪，改 SNR/出图不用重跑 GPU 仿真(~92s)。
% 若改了 STL 或 d_h_card 等仿真参数，需临时设 false 重仿真。
reuse_Mix = true;
mix_mat   = fullfile(thisdir,'put_vmd_Mix_40s.mat');

%% ==================== === 生理模型参数（与 wanzheng.m 一致） ====================
% ---------------- 呼吸（胸腔）四段式 ----------------
T_in       = 0.8;      % 吸气上升阶段 [s]
T_plateau  = 0.4;      % 平台期 [s]
T_out      = 0.8;      % 呼气下降阶段 [s]
T_b        = T_in + T_plateau + T_out + T_plateau;   % 完整呼吸周期 2.4 s
T_h_chest  = T_b;                                    % 胸腔运动周期=呼吸周期
tau_h_chest    = T_h_chest / 3;
sigma_h_chest  = T_h_chest / 5;
d_h_chest      = 0.003;      % 胸腔最大形变位移 [m]
direction_chest = [-0.5, -1.0, 0.0];

% 呼吸波形模型（put_vmd 信号级孪生）：
%   'sine'（默认）：单频平滑呼吸 0.5-0.5cos(2pi t/T_b)，周期 2.4s、零谐波，
%           使心搏带(0.8-2Hz)干净、VMD 能稳定分离心搏(1.25Hz)。雷达生命体征文献常用简化。
%   '4seg'：四段式(与 wanzheng.m 动画同)，但其高次谐波(0.83/1.25/1.67Hz…)会落在心搏带、
%           干扰心搏分离；可用 breath_smooth_s 圆滑拐角缓解。
breath_model = 'sine';
breath_smooth_s = 0.12;   % 仅 '4seg' 时使用：平滑半窗 [s]（0=不平滑）

x_min_chest = -0.13; x_max_chest = 0.13;
y_min_chest = -0.2;  y_max_chest = 0;
z_min_chest = -0.2;  z_max_chest = 0.2;

% ---------------- 心搏（心脏）高斯脉冲 ----------------
T_h_heart    = 0.8;       % 心跳周期 [s] (约 75 bpm)
tau_h_heart  = T_h_heart / 3;
sigma_h_heart= T_h_heart / 5;
% 心搏驱动模型（put_vmd 信号级孪生）：
%   'gauss'（默认）：wanzheng.m 式高斯脉冲(τ=T/3,σ=T/5)，贴"心搏搏动"；
%           在较高 SNR(如 50dB)下 VMD 可把 1.25Hz 基波干净分离成独立模态(2.5Hz 谐波在别的模态)。
%   'sine'：单频 1.25 Hz（零谐波），但会改变回波最强 bin 选择、呼吸相位谐波可能更大，较少用。
heart_model = 'gauss';
d_h_heart    = 0.003;     % 心脏最大形变位移 [m]（gauss 时：深部心脏随组织衰减见 w_heart）
direction_heart = [-0.5, -1.0, 0.0];
x_min_heart = -0.016; x_max_heart = 0.03;
y_min_heart = -0.02;  y_max_heart = 0.01;
z_min_heart = -0.022; z_max_heart = 0.0075;

use_vertex_normals_chest  = true;
use_vertex_normals_heart  = true;

% 胸腔 / 心脏反射能量权重（深层心脏经组织传播显著衰减）
w_chest = 1;
w_heart = 1e-3;

% 体表心搏微动（心前区胸壁脉冲）—— put_vmd 在胸腔表面叠加与心搏同相位的高斯形变。
% 理由：FMCW 实际测量对象是"体表复合微动"，心搏经胸壁耦合到体表（文献量级约 0.1-1 mm，
% 这里取 1 mm 以保证 VMD 能单独成模、且低 SNR 下仍稳定分离）；只靠深部心脏(权重 1e-3)
% 时心搏分量过弱无法分离。注：本改动仅作用于 put_vmd，wanzheng.m 动画不变。
use_cardiac_skin_pulse = true;
d_h_card   = 0.001;                  % 心前区体表心搏位移幅度 [m]（1 mm，文献心尖搏动上限量级；
                                     % 0.5mm 太弱 VMD 难单独成模，20dB 亦不可检出）
card_wxyz  = [0.035 0.05 0.04];      % 心前区高斯包络展宽 [m]（x, y, z）

% ---- 接收机热噪声 / 器件噪声：复高斯白噪声 n~CN(0, sigma_n^2) ----
% 在"最终回波基带"（无噪 Mix）上直接加入，即数字孪生在信号生成阶段就含噪声：
%   sigma_n = RMS(无噪回波) / 10^(SNR_dB/20)，RMS 取整个回波复数幅度均方根。
% 固定 rng 种子保证可复现；只调 SNR_dB 即可扫不同信噪比（配合 reuse_Mix=true 免重仿真）。
add_noise = true;
% 默认 50 dB：噪声按"整回波 RMS"计(σ_n=RMS/10^(SNR/20))，已尽量小但不为0（仍体现噪声项）。
% 选择说明：该 σ 基准较严苛，1mm 高斯心搏需 ~50dB 才能让 VMD 把 1.25Hz 基波干净分离成独立
% 模态（低 SNR 时 2.5Hz 谐波/噪声会混进心搏模态，频谱出现比 1.25Hz 更高的峰）。想演示低
% SNR 心搏退化，把 SNR_dB 调小即可（呼吸仍稳）。
SNR_dB    = 50;                      % 信噪比 [dB]
rng_seed  = 0;

% ---------------- 真值（由周期推出，供对照，不硬编码） ----------------
BR_true = 60 / T_b;        % ≈ 25 次/min
HR_true = 60 / T_h_heart;  % ≈ 75 次/min

%% ==================== === 雷达参数（与 wanzheng.m 一致） ====================
fc = 77.00e9;            % 载频 [Hz]
c  = 3e8;                % 光速 [m/s]
B  = 3.99034e9;          % 带宽 [Hz]
Nchirp = 2;              % 每帧 chirp 数
Nsample = 64;            % 每 chirp 采样点数
slope = 70.006e12;       % 调频斜率 [Hz/s]
Tchirp = B / slope;
sample_rate = 8000e3;    % 采样率 [Hz]
chirp_time  = Nsample/sample_rate;
idle_time   = 5e-6;
maxR  = (sample_rate*c)/(2*slope);   % 最大探测距离 [m]
value_B = chirp_time*slope;
rangeRes = c/(2*value_B);
lambda = c / fc;
dx = lambda / 2;         % RX 间距
N_rx = 8;

p_radar = [0, -5, 0];    % 雷达位置 [m]
theta_0 = 0;
phi_0   = 0;
epsilon = 1e-12;

chest_stl = fullfile(thisdir,'d.stl');   % 胸腔 STL
heart_stl = fullfile(thisdir,'b.stl');   % 心脏 STL

%% ==================== 回波仿真（可复用已保存 Mix 跳过） ====================
if reuse_Mix && exist(mix_mat,'file')
    fprintf('加载已保存 Mix: %s\n', mix_mat);
    s = load(mix_mat); Mix = s.Mix;
else
    %% ---- GPU 初始化 ----
    if gpuDeviceCount > 0
        gpu = gpuDevice(1);
        disp(['使用GPU加速: ' char(gpu.Name)]);
    else
        error('未检测到兼容GPU。请检查硬件和MATLAB Parallel Computing Toolbox。');
    end

    %% ---- 读取两个 STL 并预处理（法向量等，与 wanzheng.m 一致） ----
    [Vertices_chest, Faces_chest] = read_stl_generic(chest_stl);
    [Vertices_heart,  Faces_heart ] = read_stl_generic(heart_stl);
    Vertices_chest = double(Vertices_chest); Faces_chest = double(Faces_chest);
    Vertices_heart = double(Vertices_heart); Faces_heart = double(Faces_heart);
    numVerts_chest = size(Vertices_chest,1);
    numVerts_heart = size(Vertices_heart,1);

    % chest vertex normals
    numFacesC = size(Faces_chest,1);
    n_sum_c = zeros(numVerts_chest,3);
    for f = 1:numFacesC
        i1 = Faces_chest(f,1); i2 = Faces_chest(f,2); i3 = Faces_chest(f,3);
        v1 = Vertices_chest(i1,:); v2 = Vertices_chest(i2,:); v3 = Vertices_chest(i3,:);
        nf = cross(v2 - v1, v3 - v1); nf = nf/(norm(nf)+epsilon);
        n_sum_c(i1,:)=n_sum_c(i1,:)+nf; n_sum_c(i2,:)=n_sum_c(i2,:)+nf; n_sum_c(i3,:)=n_sum_c(i3,:)+nf;
    end
    vertex_normals_chest = n_sum_c ./ (sqrt(sum(n_sum_c.^2,2))+epsilon);

    % heart vertex normals
    numFacesH = size(Faces_heart,1);
    n_sum_h = zeros(numVerts_heart,3);
    for f = 1:numFacesH
        i1 = Faces_heart(f,1); i2 = Faces_heart(f,2); i3 = Faces_heart(f,3);
        v1 = Vertices_heart(i1,:); v2 = Vertices_heart(i2,:); v3 = Vertices_heart(i3,:);
        nf = cross(v2 - v1, v3 - v1); nf = nf/(norm(nf)+epsilon);
        n_sum_h(i1,:)=n_sum_h(i1,:)+nf; n_sum_h(i2,:)=n_sum_h(i2,:)+nf; n_sum_h(i3,:)=n_sum_h(i3,:)+nf;
    end
    vertex_normals_heart = n_sum_h ./ (sqrt(sum(n_sum_h.^2,2))+epsilon);

    %% ---- 形变幅度 A_i 与方向 delta_max（与 wanzheng.m 一致） ----
    e_const = exp(1);
    A_chest = zeros(numVerts_chest,1);
    for i=1:numVerts_chest
        x=Vertices_chest(i,1); y=Vertices_chest(i,2); z=Vertices_chest(i,3);
        if ~(x>=x_min_chest && x<=x_max_chest && y>=y_min_chest && y<=y_max_chest && z>=z_min_chest && z<=z_max_chest)
            continue;
        end
        fz = -( (z-z_min_chest)*(z-z_max_chest) )/( (z_max_chest-z_min_chest)^2 ); fz=max(fz,0);
        fy = -( (y-y_min_chest)*(y-y_max_chest) )/( (y_max_chest-y_min_chest)^2 ); fy=max(fy,0);
        sqrt_fz = fz;                 % 胸腔原代码这里无 sqrt
        sqrt_fy = sqrt(fy);
        ratio = (x-x_min_chest)/(x_max_chest-x_min_chest); ratio=min(max(ratio,0),1);
        log_term = log( e_const + (1-e_const)*ratio );
        A_chest(i) = d_h_chest * sqrt_fz * 1 * log_term;
    end

    A_heart = zeros(numVerts_heart,1);
    for i=1:numVerts_heart
        x=Vertices_heart(i,1); y=Vertices_heart(i,2); z=Vertices_heart(i,3);
        if ~(x>=x_min_heart && x<=x_max_heart && y>=y_min_heart && y<=y_max_heart && z>=z_min_heart && z<=z_max_heart)
            continue;
        end
        fz = -( (z-z_min_heart)*(z-z_max_heart) )/( (z_max_heart-z_min_heart)^2 ); fz=max(fz,0);
        fy = -( (y-y_min_heart)*(y-y_max_heart) )/( (y_max_heart-y_min_heart)^2 ); fy=max(fy,0);
        sqrt_fz = sqrt(fz);
        sqrt_fy = sqrt(fy);
        ratio = (x-x_min_heart)/(x_max_heart-x_min_heart); ratio=min(max(ratio,0),1);
        log_term = log( e_const + (1-e_const)*ratio );
        A_heart(i) = d_h_heart * sqrt_fz * sqrt_fy * log_term;
    end

    % 方向
    if use_vertex_normals_chest
        delta_max_chest = A_chest .* vertex_normals_chest;
    else
        dir_c = direction_chest/norm(direction_chest);
        delta_max_chest = A_chest .* dir_c;
    end
    if use_vertex_normals_heart
        delta_max_heart = A_heart .* vertex_normals_heart;
    else
        dir_h = direction_heart/norm(direction_heart);
        delta_max_heart = A_heart .* dir_h;
    end

    % 合并顶点与两个形变场（GPU 上按"胸腔块 + 心脏块"统一向量计算）
    %   breathField：呼吸场，仅胸腔表面（随 M_chest 缩放）
    %   cardField  ：心搏场 = 心前区胸壁脉冲(d_h_card) + 深部心脏(d_max_heart)，随 M_heart 缩放
    [~, pk] = max(A_chest);                      % 心前区锚点：胸腔形变峰值顶点
    cax = Vertices_chest(pk,1); cay = Vertices_chest(pk,2); caz = Vertices_chest(pk,3);
    if use_cardiac_skin_pulse
        g_card = zeros(numVerts_chest,1);
        for i = 1:numVerts_chest
            if A_chest(i) <= 0, continue; end     % 仅在胸腔形变区（前胸壁）加脉冲
            dxc = (Vertices_chest(i,1)-cax)/card_wxyz(1);
            dyc = (Vertices_chest(i,2)-cay)/card_wxyz(2);
            dzc = (Vertices_chest(i,3)-caz)/card_wxyz(3);
            g_card(i) = exp(-(dxc^2 + dyc^2 + dzc^2));
        end
        delta_card_chest = d_h_card * g_card .* vertex_normals_chest;
    else
        delta_card_chest = zeros(numVerts_chest,3);
    end

    Vertices_all = [Vertices_chest; Vertices_heart];
    vertex_normals_all = [vertex_normals_chest; vertex_normals_heart];
    breathField = [delta_max_chest; zeros(numVerts_heart,3)];
    cardField   = [delta_card_chest; delta_max_heart];

    Vertices_all_gpu = gpuArray(Vertices_all);
    breathF_gpu = gpuArray(breathField);
    cardF_gpu   = gpuArray(cardField);
    vertex_normals_all_gpu = gpuArray(vertex_normals_all);
    V_prev_gpu = Vertices_all_gpu;

    t_frame_fmcw_gpu = gpuArray(linspace(0, Nchirp*Tchirp, Nsample*Nchirp));

    Mix = zeros(N_rx, Nsample, Nchirp, Nframe);
    fprintf('开始 %d s / %d 帧 回波仿真（进度每 50 帧打印一次）...\n', total_seconds, Nframe);
    tSim = tic;
    % 呼吸单周期波形（按帧查表，节律仍 T_b=2.4 s = 25 次/min）
    perB = round(T_b * fps);                        % 每周期帧数（=60）
    ttB  = (0:perB-1)/fps;
    switch lower(breath_model)
        case 'sine'
            MpB = 0.5*(1 - cos(2*pi*ttB/T_b));      % 单频、零谐波
        case '4seg'
            MpB = zeros(1, perB);
            for j = 1:perB
                tmc = ttB(j);
                if tmc < T_in
                    MpB(j) = sin( (pi/2)*(tmc/T_in) )^2;
                elseif tmc < T_in + T_plateau
                    MpB(j) = 1;
                elseif tmc < T_in + T_plateau + T_out
                    xq = (tmc - T_in - T_plateau)/T_out;
                    MpB(j) = sin( (pi/2) + (pi/2)*xq )^2;
                else
                    MpB(j) = 0;
                end
            end
            if breath_smooth_s > 0                  % 可选：圆周平滑削拐角谐波
                nker = max(1, round(breath_smooth_s*fps));
                kk   = -nker:nker;
                ker  = (1 + cos(pi*kk/nker))/2; ker = ker/sum(ker);
                MpB = cconv(MpB, ker, perB);
            end
        otherwise
            error('breath_model 需为 ''sine'' 或 ''4seg''');
    end
    for frame = 1:Nframe
        M_t_chest = MpB(1 + mod(frame-1, perB));    % 呼吸四段式（平滑后）
        t = (frame-1)/fps;
        % 心搏驱动（心脏/体表心搏）
        t_mod_heart = mod(t, T_h_heart);
        switch lower(heart_model)
            case 'sine'
                M_t_heart = sin(2*pi*t_mod_heart/T_h_heart);       % 单频 1.25 Hz
            case 'gauss'
                M_t_heart = exp(-((t_mod_heart - tau_h_heart)/sigma_h_heart)^2);
            otherwise
                error('heart_model 需为 ''sine'' 或 ''gauss''');
        end

        V_frame_gpu = Vertices_all_gpu + M_t_chest*breathF_gpu + M_t_heart*cardF_gpu;

        % 顶点速度（帧间位移 / dt）
        v_obj_gpu = (V_frame_gpu - V_prev_gpu) * fps;   % fps = 1/dt
        V_prev_gpu = V_frame_gpu;

        % FMCW 回波合成（逐 RX 累加）
        for i_rx = 1:N_rx
            p_rx = gpuArray([p_radar(1), p_radar(2)+(i_rx-1)*dx, p_radar(3)]);
            d_i_gpu = V_frame_gpu - p_rx;
            d_i_norm_gpu = sqrt(sum(d_i_gpu.^2,2));
            d_hat_gpu = d_i_gpu ./ (d_i_norm_gpu + epsilon);
            v_i_gpu = sum(d_i_gpu .* v_obj_gpu,2) ./ (d_i_norm_gpu + epsilon);
            s_i_gpu = max(0, -sum(d_hat_gpu .* vertex_normals_all_gpu,2)); % Lambertian
            s_weights = [w_chest*ones(numVerts_chest,1,'gpuArray'); ...
                         w_heart*ones(numVerts_heart,1,'gpuArray')];
            s_i_gpu = s_i_gpu .* s_weights;

            r_t_gpu = d_i_norm_gpu + v_i_gpu .* t_frame_fmcw_gpu;  % 广播
            td_gpu  = 2 * r_t_gpu / c;
            phi_tx_gpu = 2*pi*( fc*t_frame_fmcw_gpu + 0.5*slope*t_frame_fmcw_gpu.^2 );
            phi_rx_gpu = 2*pi*( fc*(t_frame_fmcw_gpu - td_gpu) + 0.5*slope*(t_frame_fmcw_gpu - td_gpu).^2 );
            mix_point_gpu = exp(1j*phi_tx_gpu) .* conj(exp(1j*phi_rx_gpu).*s_i_gpu);
            mix_frame = gather(sum(mix_point_gpu,1));
            Mix(i_rx,:,:,frame) = reshape(mix_frame, Nsample, Nchirp);
        end

        if mod(frame,50)==0 || frame==1
            fprintf('  frame %4d/%d  (M_chest=%.3f M_heart=%.3f)  t=%.1fs\n', ...
                    frame, Nframe, M_t_chest, M_t_heart, toc(tSim));
        end
    end
    fprintf('回波仿真完成: %.1f s（%d 帧 x %d chirp x %d 采样）\n', ...
            toc(tSim), Nframe, Nchirp, Nsample);

    if ~reuse_Mix
        save(mix_mat,'Mix','-v7.3');
        fprintf('Mix 已保存至 %s（下次可设 reuse_Mix=true 跳过仿真）\n', mix_mat);
    end
end

[N_rx, Nsample, Nchirp, Nframe] = size(Mix);   % 复用 Mix 时也能取对维度

%% ==================== 接收机噪声：最终回波 = 无噪 Mix + 复高斯白噪声 ====================
% 现实雷达系统受热噪声与器件噪声影响：n(t)~CN(0,sigma_n^2)，
%   sigma_n = RMS(R^(c)(tau,t)) / 10^(SNR/20)，
%   R~ = R + n(t)。这里 R 即仿真得到的无噪基带回波 Mix（数字孪生生成阶段即含噪）。
if add_noise
    rng(rng_seed);                                   % 保证可复现
    R_rms   = sqrt(mean(abs(Mix(:)).^2));            % 无噪回波复数幅度 RMS
    sigma_n = R_rms / 10^(SNR_dB/20);
    Mix = Mix + (sigma_n/sqrt(2)) * (randn(size(Mix)) + 1i*randn(size(Mix)));
    fprintf('已加入复高斯白噪声: SNR=%.1f dB, sigma_n=%.3e (RMS/10^(SNR/20))\n', SNR_dB, sigma_n);
else
    fprintf('未加噪声（add_noise=false）。\n');
end

%% ==================== 距离维 FFT 与慢时间相位提取 ====================
% 取 RX1 第一 chirp（与 main_test / wanzheng 注释段一致）
mix1 = squeeze(Mix(1,:,1,:));                 % Nsample x Nframe
RangeFFT = fft(mix1 .* hanning(Nsample), Nsample, 1);   % 加汉宁窗后沿快时间做 FFT

% 自动选最强距离 bin（主目标所在 bin）
[~, MaxIndex] = max(mean(abs(RangeFFT),2));
fprintf('选中距离 bin %d（≈ %.3f m，主目标胸腔）\n', MaxIndex, (MaxIndex-1)*rangeRes);

% 相位提取：直接对该 bin 的复数相量 angle 做 unwrap 后逐帧差分（不去慢时间均值）。
% 去均值会把相量中心移到原点附近，噪声下增量相位偶发尖刺污染低频谱（因此 main_test
% 实测里才要 f_phaseDenoise 去刺）；仿真主目标相量幅值恒定很大，不靠近原点，
% 直接 unwrap 对噪声更稳健，呼吸率在低 SNR 下仍稳定。
phase_bin = unwrap(angle(RangeFFT(MaxIndex,:)));      % 慢时间相位 [rad]
angleDenoised = diff(phase_bin);                      % 相位增量序列（Nframe-1 点）

% 备选：main_test 的去均值 + 增量相位公式（保留供对照）
% RangeFFTcor = RangeFFT - mean(RangeFFT,2);
% I = real(RangeFFTcor(MaxIndex,:)); Q = imag(RangeFFTcor(MaxIndex,:));
% phi = zeros(1,Nframe);
% for k=2:Nframe
%     dI=I(k)-I(k-1); dQ=Q(k)-Q(k-1);
%     phi(k)=phi(k-1)+(I(k)*dQ-Q(k)*dI)/(I(k)^2+Q(k)^2);
% end
% angleDenoised = diff(phi);

% 实测数据才做的脉冲去噪（main_test 用 f_phaseDenoise，仿真数据无脉冲噪声故省略）
% angleDenoised2 = angleDenoised;
% for ii=3:numel(angleDenoised)
%     angleDenoised2(ii)=f_phaseDenoise(angleDenoised(ii-2),angleDenoised(ii-1),angleDenoised(ii),0.3);
% end

%% ==================== VMD 变分模态分解（参数同 main_test.m） ====================
K     = 7;       % 模态数
alpha = 1000;    % 惩罚因子
tol   = 1e-6;    % 收敛容差
[imf, IMFs_fft_raw, ~] = VMD(angleDenoised, alpha, 0, K, 0, 2, tol);

IMFs_fft = abs(IMFs_fft_raw)';                % K x M（每行一个模态频谱，居中）
M        = size(IMFs_fft,2);
freqs    = (-M/2:M/2-1)*(fs/M);               % 居中频率轴 [Hz]

% 图注/输出通用时间轴（VMD 输出长度可能与原信号一致或补一位）
t_axis = (0:size(imf,2)-1)/fs;

%% ==================== 自动挑选 呼吸 / 心搏 IMF（频段同 main_test） ====================
resp_range = [0.1 0.5];   % 呼吸频段 [Hz]
hb_range   = [0.8 2.0];   % 心搏频段 [Hz]

resp_amp = zeros(1,K); hb_amp = zeros(1,K);
for i = 1:K
    resp_amp(i) = max(IMFs_fft(i, freqs>=resp_range(1) & freqs<=resp_range(2)));
    hb_amp(i)   = max(IMFs_fft(i, freqs>=hb_range(1)   & freqs<=hb_range(2)));
end
[~, resp_idx] = max(resp_amp);      % 呼吸主 IMF
[~, hb_idx]   = max(hb_amp);        % 心搏主 IMF

resp_IMF  = imf(resp_idx,:);        % 呼吸 IMF（时域）
hb_IMF    = imf(hb_idx,:);          % 心搏 IMF（时域）
breath_fre = IMFs_fft(resp_idx,:);  % 呼吸 IMF 频谱（居中）
heart_fre  = IMFs_fft(hb_idx,:);    % 心搏 IMF 频谱（居中）

fprintf('VMD 选出：呼吸主IMF = imf%d，心搏主IMF = imf%d\n', resp_idx, hb_idx);

%% ==================== 呼吸率 / 心率估计（谱峰->次/分） ====================
maskR = freqs>=resp_range(1) & freqs<=resp_range(2);
[~,ir] = max(breath_fre .* maskR);    BR_meas = freqs(ir)*60;   % bpm
maskH = freqs>=hb_range(1) & freqs<=hb_range(2);
[~,ih] = max(heart_fre  .* maskH);    HR_meas = freqs(ih)*60;   % bpm

% 交叉验证：心搏用 pwelch 谱估计（同 main_test）
[Pxx_h, fh] = pwelch(hb_IMF(:)-mean(hb_IMF), [], [], [], fs);
maskpw = fh>=hb_range(1) & fh<=hb_range(2);
[~,ipw] = max(Pxx_h(maskpw));
HR_pw = fh(maskpw); HR_pw = HR_pw(ipw)*60;

fprintf('\n========== VMD 结果 ==========\n');
fprintf('呼吸率 BR: 测得 %.2f 次/min   (真值 %.2f 次/min)\n', BR_meas, BR_true);
fprintf('心率   HR: 测得 %.2f 次/min   (真值 %.2f 次/min)\n', HR_meas, HR_true);
fprintf('心搏 pwelch 交叉验证 HR = %.2f 次/min\n', HR_pw);
fprintf('==============================\n');

%% ==================== 出图：全 IMF（时域） ====================
fig = figure('Color','w','Position',[40 40 850 1200]);
for k = 1:K
    subplot(K,1,k);
    plot(t_axis, imf(k,:),'LineWidth',0.8); grid on;
    title(['IMF',num2str(k)]);
    ylabel('幅度');
    xlim([0 t_axis(end)]);
    if k==K, xlabel('时间 (s)'); end
end
sgtitle('VMD 分解的 7 个本征模态（时域）');
exportgraphics(fig, fullfile(outdir,'fig6_VMD分解_全IMF时域.png'),'Resolution',300);
close(fig);

%% ==================== 出图：全 IMF（频谱 0-5 Hz） ====================
fig = figure('Color','w','Position',[40 40 850 1200]);
for k = 1:K
    subplot(K,1,k);
    plot(freqs, IMFs_fft(k,:),'LineWidth',0.8); grid on;
    title(['IMF',num2str(k),' 频谱']);
    xlim([0 5]); ylabel('幅度');
    if k==K, xlabel('频率 (Hz)'); end
end
sgtitle('VMD 各本征模态频谱（0–5 Hz）');
exportgraphics(fig, fullfile(outdir,'fig7_VMD分解_全IMF频谱.png'),'Resolution',300);
close(fig);

%% ==================== 出图：呼吸/心搏分离（2x2，同 main_test VMD结果图） ====================
[~,irp] = max(breath_fre .* maskR); fr_peak = freqs(irp);
[~,ihp] = max(heart_fre  .* maskH); fh_peak = freqs(ihp);

fig = figure('Color','w','Position',[80 80 1050 800]);
subplot(2,2,1);
plot(t_axis, resp_IMF,'LineWidth',0.9,'Color',[0 0.45 0.75]); grid on;
title(['呼吸信号时域（IMF',num2str(resp_idx),'）']); ylabel('幅度'); xlabel('时间 (s)');
subplot(2,2,2);
plot(freqs, breath_fre,'LineWidth',0.9,'Color',[0 0.45 0.75]); grid on; hold on;
plot(fr_peak, max(breath_fre(maskR)),'rv','MarkerSize',8,'LineWidth',1.5);
title(['呼吸信号频谱（IMF',num2str(resp_idx),'）']); ylabel('幅度'); xlabel('频率 (Hz)'); xlim([0 3]);
text(fr_peak, max(breath_fre(maskR)), sprintf('  %.2f Hz = %.1f 次/min',fr_peak,BR_meas),'Color','r');
subplot(2,2,3);
plot(t_axis, hb_IMF,'LineWidth',0.9,'Color',[0.8 0.2 0.2]); grid on;
title(['心跳信号时域（IMF',num2str(hb_idx),'）']); ylabel('幅度'); xlabel('时间 (s)');
subplot(2,2,4);
plot(freqs, heart_fre,'LineWidth',0.9,'Color',[0.8 0.2 0.2]); grid on; hold on;
plot(fh_peak, max(heart_fre(maskH)),'rv','MarkerSize',8,'LineWidth',1.5);
title(['心跳信号频谱（IMF',num2str(hb_idx),'）']); ylabel('幅度'); xlabel('频率 (Hz)'); xlim([0 2]);
text(fh_peak, max(heart_fre(maskH)), sprintf('  %.2f Hz = %.1f 次/min',fh_peak,HR_meas),'Color','r');
sgtitle('VMD 分离的呼吸与心搏信号（时域 + 频谱）');
exportgraphics(fig, fullfile(outdir,'fig8_VMD呼吸心搏分离_时域频谱.png'),'Resolution',300);
close(fig);

%% ==================== 出图：测得/真值对照（总体系一归一化，看实际相对幅度） ====================
% 不各自归一化：呼吸、心搏共用同一个归一化尺度 max(全谱)，保留两者真实幅度比。
common_scale = max([breath_fre, heart_fre]);
[~, ir9] = max(breath_fre .* maskR);
[~, ih9] = max(heart_fre  .* maskH);
amp_ratio_hb2br = heart_fre(ih9) / breath_fre(ir9);   % 心搏峰幅度 / 呼吸峰幅度

fig = figure('Color','w','Position',[100 100 950 600]);
plot(freqs, breath_fre/common_scale,'b','LineWidth',1.2); hold on;
plot(freqs, heart_fre/common_scale,'r','LineWidth',1.2);
plot(freqs(ir9), breath_fre(ir9)/common_scale,'bo','MarkerFaceColor','b');
plot(freqs(ih9), heart_fre(ih9)/common_scale,'ro','MarkerFaceColor','r');
xlim([0 3]); grid on;
xlabel('频率 (Hz)'); ylabel('归一化幅度（总体系一，保留实际幅度比）');
legend({'呼吸 IMF','心搏 IMF'},'Location','northeast');
title('VMD 分离模态谱峰对比（总体归一化，未各自归一化）');
dim = [0.57 0.14 0.42 0.24];
str = sprintf(['呼吸率 BR：测得 %.1f 次/min  真值 %.1f\n', ...
               '心率   HR：测得 %.1f 次/min  真值 %.1f\n', ...
               '谱峰幅度比 心搏/呼吸 = %.3f'], ...
              BR_meas, BR_true, HR_meas, HR_true, amp_ratio_hb2br);
annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','w','EdgeColor','k');
exportgraphics(fig, fullfile(outdir,'fig9_VMD呼吸心率估计_真值对照.png'),'Resolution',300);
close(fig);

%% ==================== 独立出图：逐模态 时域 / 频谱（每模态各一张） ====================
for k = 1:K
    fig = figure('Color','w','Position',[40 40 820 360]);
    plot(t_axis, imf(k,:),'LineWidth',0.8); grid on;
    title(['IMF',num2str(k),' 时域']); xlabel('时间 (s)'); ylabel('幅度');
    exportgraphics(fig, fullfile(outdir,['vmd_imf',num2str(k),'_时域.png']),'Resolution',300);
    close(fig);

    fig = figure('Color','w','Position',[40 40 820 360]);
    plot(freqs, IMFs_fft(k,:),'LineWidth',0.8); grid on;
    title(['IMF',num2str(k),' 频谱']); xlabel('频率 (Hz)'); ylabel('幅度'); xlim([0 5]);
    exportgraphics(fig, fullfile(outdir,['vmd_imf',num2str(k),'_频谱.png']),'Resolution',300);
    close(fig);
end

%% ==================== 独立出图：呼吸/心搏 时域与频谱（各一张，即 fig8 拆开） ====================
% 呼吸 时域
fig = figure('Color','w','Position',[40 40 820 360]);
plot(t_axis, resp_IMF,'LineWidth',0.9,'Color',[0 0.45 0.75]); grid on;
title(['呼吸信号时域（IMF',num2str(resp_idx),'）']); xlabel('时间 (s)'); ylabel('幅度');
exportgraphics(fig, fullfile(outdir,'vmd_呼吸_时域.png'),'Resolution',300);
close(fig);
% 呼吸 频谱（含谱峰标注）
fig = figure('Color','w','Position',[40 40 820 360]);
plot(freqs, breath_fre,'LineWidth',0.9,'Color',[0 0.45 0.75]); grid on; hold on;
plot(freqs(ir9), breath_fre(ir9),'rv','MarkerSize',8,'LineWidth',1.5);
title(['呼吸信号频谱（IMF',num2str(resp_idx),'）']); xlabel('频率 (Hz)'); ylabel('幅度'); xlim([0 3]);
text(freqs(ir9), breath_fre(ir9), sprintf('  %.2f Hz = %.1f 次/min',freqs(ir9),BR_meas),'Color','r');
exportgraphics(fig, fullfile(outdir,'vmd_呼吸_频谱.png'),'Resolution',300);
close(fig);% 心跳 时域
fig = figure('Color','w','Position',[40 40 820 360]);
plot(t_axis, hb_IMF,'LineWidth',0.9,'Color',[0.8 0.2 0.2]); grid on;
title(['心跳信号时域（IMF',num2str(hb_idx),'）']); xlabel('时间 (s)'); ylabel('幅度');
exportgraphics(fig, fullfile(outdir,'vmd_心跳_时域.png'),'Resolution',300);
close(fig);
% 心跳 频谱（含谱峰标注）
fig = figure('Color','w','Position',[40 40 820 360]);
plot(freqs, heart_fre,'LineWidth',0.9,'Color',[0.8 0.2 0.2]); grid on; hold on;
plot(freqs(ih9), heart_fre(ih9),'rv','MarkerSize',8,'LineWidth',1.5);
title(['心跳信号频谱（IMF',num2str(hb_idx),'）']); xlabel('频率 (Hz)'); ylabel('幅度'); xlim([0 2]);
text(freqs(ih9), heart_fre(ih9), sprintf('  %.2f Hz = %.1f 次/min',freqs(ih9),HR_meas),'Color','r');
exportgraphics(fig, fullfile(outdir,'vmd_心跳_频谱.png'),'Resolution',300);
close(fig);

%% ==================== 独立出图：雷达 距离-时间图（RX1 第一 chirp，加噪后回波） ====================
range_axis = (0:Nsample-1) * rangeRes;               % 距离 [m]
time_axis2 = (0:Nframe-1) / fs;                       % 时间 [s]
fig = figure('Color','w','Position',[100 100 900 520]);
mesh(range_axis, time_axis2, abs(RangeFFT).'/max(abs(RangeFFT(:))));
view([0 90]); colorbar;
xlabel('距离 (m)'); ylabel('时间 (s)');
title(['雷达距离-时间图（RX1 第一chirp，加噪 SNR=',num2str(SNR_dB),' dB）']);
exportgraphics(fig, fullfile(outdir,'fig10_雷达距离时间图.png'),'Resolution',300);
close(fig);

fprintf('图片已保存到：\n  %s\n', outdir);
fprintf('全部完成。\n');

%% ==================== 本地函数：通用 STL 读取（与 wanzheng.m 相同） ====================
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
        V = unique(rawV,'rows','stable');
        F = rawF;
    else
        fseek(fid,80,'bof');
        numFaces = fread(fid,1,'uint32');
        rawV = zeros(numFaces*3,3);
        rawF = reshape(1:(numFaces*3),3,[])';
        for k=1:numFaces
            fread(fid,3,'float32');           % normal
            v1 = fread(fid,3,'float32')';
            v2 = fread(fid,3,'float32')';
            v3 = fread(fid,3,'float32')';
            rawV((k-1)*3+1,:) = v1;
            rawV((k-1)*3+2,:) = v2;
            rawV((k-1)*3+3,:) = v3;
            fread(fid,1,'uint16');            % attr
        end
        fclose(fid);
        [V,~,idxMap] = uniquetol(rawV, 1e-9, 'ByRows', true);
        F = reshape(idxMap, 3, [])';
    end
end
