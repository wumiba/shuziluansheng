clear; close all; clc;

%% ---------------------- 1. 参数设置 ----------------------
% 呼吸物理参数
T_in = 0.8; T_plateau = 0.4; T_out = 0.8; 
T_b = T_in + T_plateau + T_out + T_plateau; 
T_h = T_b; 
fps = 100;            % 高帧率保证 unwrap 成功
d_h = 0.08;           % 最大形变位移 [m]
total_seconds = 8;
total_frames = round(total_seconds * fps);

% 心脏/呼吸影响范围（空间形变选通框）
daxiao = 1;
x_min = -0.13*daxiao; x_max = 0.13*daxiao;
y_min = -0.2*daxiao;  y_max = 0;
z_min = -0.2*daxiao;  z_max = 0.2*daxiao;

% 雷达 FMCW 参数
fc = 77.00e9;             
c = 3e8;                  
lambda = c / fc;
B = 3.99e9;            
Nchirp = 1;               
Nsample = 128;            
slope = 70.006e12;        
Tchirp = B / slope;       
sample_rate = 8000e3;     
maxR = (sample_rate*c)/(2*slope);    

% 雷达位置
p_radar = [0, -0.5, 0];   
epsilon = 1e-12;

%% ---------------------- 2. 初始化 GPU 与 模型 ----------------------
if gpuDeviceCount > 0
    gpu = gpuDevice(1);
    disp('✅ GPU 加速已开启: ' + string(gpu.Name));
else
    error('❌ 需要 GPU 支持以运行高密度散射点仿真');
end

% STL 文件读取（请确保当前目录有 d.stl）
stl_file = 'd.stl'; 
[Vertices, Faces] = read_stl_generic(stl_file);
Vertices = double(Vertices);
numVerts = size(Vertices,1);

% 计算顶点法向量
n_sum = zeros(numVerts,3);
for f = 1:size(Faces,1)
    idx = Faces(f,:);
    v1 = Vertices(idx(1),:); v2 = Vertices(idx(2),:); v3 = Vertices(idx(3),:);
    nf = cross(v2 - v1, v3 - v1);
    nf = nf / (norm(nf) + epsilon);
    n_sum(idx,:) = n_sum(idx,:) + nf;
end
vertex_normals = n_sum ./ (sqrt(sum(n_sum.^2,2)) + epsilon);

% 计算形变权重 A (仅在选通框内的点会动)
A = zeros(numVerts,1);
for i=1:numVerts
    v = Vertices(i,:);
    if (v(1)>=x_min && v(1)<=x_max && v(2)>=y_min && v(2)<=y_max && v(3)>=z_min && v(3)<=z_max)
        A(i) = d_h * exp(-((v(1)^2 + v(3)^2)/0.01)); 
    end
end

% 将基础数据移入 GPU
V_base_gpu = gpuArray(Vertices);
N_base_gpu = gpuArray(vertex_normals);
A_gpu = gpuArray(A);
p_radar_gpu = gpuArray(p_radar);
t_fast_gpu = gpuArray(linspace(0, Tchirp, Nsample)); 

%% ---------------------- 3. 动画与信号仿真循环 (无空间选通) ----------------------
Mix = zeros(Nsample, total_frames); 
V_prev_gpu = V_base_gpu;

% 预分配绘图
figure('Color','w','Position',[100 100 1000 400]);
subplot(1,2,1);
h_patch = patch('Vertices', Vertices, 'Faces', Faces, 'FaceColor', [0.8 0.3 0.3], 'EdgeColor', 'none');
axis equal; view(3); camlight; lighting gouraud; grid on;
title('3D 形变仿真 (全顶点散射)');

for frame = 1:total_frames
    t = (frame - 1) / fps;
    t_mod = mod(t, T_h);
    
    % --- 呼吸调制函数 M_t ---
    if t_mod < T_in
        M_t = sin((pi/2)*(t_mod/T_in))^2;
    elseif t_mod < T_in + T_plateau
        M_t = 1;
    elseif t_mod < T_in + T_plateau + T_out
        M_t = sin((pi/2) + (pi/2)*((t_mod-T_in-T_plateau)/T_out))^2;
    else
        M_t = 0;
    end
    
    % --- 更新顶点位置与速度 ---
    V_curr_gpu = V_base_gpu + (M_t .* A_gpu .* N_base_gpu);
    v_obj_gpu = (V_curr_gpu - V_prev_gpu) * fps;
    V_prev_gpu = V_curr_gpu;
    
    % --- 雷达几何计算 (无 Mask, 计算所有点) ---
    vec_to_radar = p_radar_gpu - V_curr_gpu; 
    dist_gpu = sqrt(sum(vec_to_radar.^2, 2));         % 距离 r0
    unit_vec = vec_to_radar ./ (dist_gpu + epsilon);    % 方向
    
    v_rad = -sum(unit_vec .* v_obj_gpu, 2);           % 径向速度
    
    % 反射增益：基于 Lambertian 模型 (法线与雷达夹角)
    cos_theta = max(0, sum(unit_vec .* N_base_gpu, 2)); 
    s_gain = cos_theta ./ (dist_gpu.^2 + epsilon); 
    
    % --- 信号合成 (GPU 广播计算) ---
    r_t = dist_gpu + v_rad * t_fast_gpu; 
    td = 2 * r_t / c;
    
    % FMCW 混频信号相位计算
    phi = 2 * pi * (fc * td + slope * t_fast_gpu .* td);
    mix_points = s_gain .* exp(1j * phi);
    
    % 累加所有散射点贡献
    Mix(:, frame) = gather(sum(mix_points, 1));
    
    % 更新 3D 动画可视化
    if mod(frame, 10) == 1
        set(h_patch, 'Vertices', gather(V_curr_gpu));
        drawnow limitrate;
    end
end

%% ---------------------- 4. 信号处理与结果展示 ----------------------
% 1. 距离 FFT
range_fft = fft(Mix, [], 1);
range_abs = abs(range_fft(1:Nsample/2, :));

% 2. 动态定位目标距离 Bin (逐帧找最大能量点)
[~, max_idx_per_frame] = max(range_abs, [], 1);
% 为了稳定，可以取全局平均能量最大的 Bin 作为参考，或者逐帧提取
[~, global_max_idx] = max(mean(range_abs, 2));

phase_series = zeros(1, total_frames);
for f = 1:total_frames
    phase_series(f) = angle(range_fft(global_max_idx, f));
end

% 3. 相位解卷绕与去趋势
unwrapped_phase = unwrap(phase_series);
detrended_phase = detrend(unwrapped_phase); 

% --- 绘图 A: 相位曲线 ---
subplot(1,2,2);
t_axis = (0:total_frames-1)/fps;
plot(t_axis, detrended_phase, 'LineWidth', 1.5, 'Color', [0 0.447 0.741]);
grid on;
xlabel('时间 (s)'); ylabel('相位 (rad)');
title(['提取的呼吸相位 (距离 Bin: ', num2str(global_max_idx), ')']);

% --- 绘图 B: RTI 能量图 ---
figure('Name','Range-Time Intensity (RTI)','Color','w');
dist_axis = (0:Nsample/2-1) * (maxR / (Nsample/2));
imagesc(t_axis, dist_axis, range_abs);
axis xy; colormap('jet'); colorbar;
xlabel('时间 (s)'); ylabel('距离 (m)');
title('RTI 能量图 (全模型散射)');

