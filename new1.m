% heartbeat_stl_animation_with_normals_and_radar.m
% 心脏 STL 形变动画 + 雷达回波信号建模（多目标、多通道、Lambertian反射 + 高斯噪声）

clear; close all; clc;

%% ---------------------- 用户参数 ----------------------
T_h      = 0.8;        % 心跳周期（秒）
fps      = 25;         % 帧率
tau_h    = T_h / 3;    % 心跳强度最大时刻
sigma_h  = T_h / 10;   % 控制心跳尖锐程度

d_h      = 0.003;      % 最大形变位移量
direction = [-0.5, -1.0, 0.0];  % 固定形变方向向量

% 心脏影响范围
x_min = -0.016; x_max = 0.03;
y_min = -0.02; y_max = 0.01;
z_min = -0.022; z_max = 0.0075;

% 动画长度
total_seconds = 15;
total_frames = round(total_seconds * fps);

% STL 文件路径
stl_file = 'b.stl';

% 是否使用顶点法向量
use_vertex_normals = true;%true以法向量形式扩张，false沿固定方向扩张



%% ---------------------- 读取 STL ----------------------
[Vertices, Faces] = read_stl_generic(stl_file);
Vertices = double(Vertices);
Faces = double(Faces);

%% ---------------------- 计算顶点法向量 ----------------------
epsilon = 1e-12;
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
    sqrt_fz = sqrt(fz);
    sqrt_fy = sqrt(fy);
    ratio = (x - x_min) / (x_max - x_min);
    ratio = min(max(ratio,0),1);
    log_term = log( e_const + (1 - e_const) * ratio );
    A(i) = d_h * sqrt_fz * sqrt_fy * log_term;
end

%% ---------------------- 确定形变方向 ----------------------
if use_vertex_normals
    delta_max = A .* vertex_normals;
else
    direction = direction / norm(direction); % 归一化
    delta_max = A .* direction;
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
xlim(ax, vc(1) + [-0.05 0.05]);
ylim(ax, vc(2) + [-0.05 0.05]);
zlim(ax, vc(3) + [-0.05 0.05]);
drawnow;

%% ---------------------- 动画 ----------------------
for frame = 0:(total_frames-1)
    t = frame / fps;
    t_mod = mod(t, T_h);
    M_t = exp(-((t_mod - tau_h) / sigma_h)^2);
    V_frame = Vertices + (M_t .* delta_max);
    set(p, 'Vertices', V_frame);
    title(ax, sprintf('Heart Apex Beating — frame %d / %d (M_t=%.3f)', ...
        frame+1, total_frames, M_t));
    drawnow;
    pause(1/fps);
end

fprintf('✅ 动画完成: %d 帧, %d FPS\n', total_frames, fps);


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
        vertexTokens = regexp(str, 'vertex\s+([-+eE0-9\.\-]+)\s+([-+eE0-9\.\-]+)\s+([-+eE0-9\.\-]+)', 'tokens');
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
