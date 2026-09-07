model = stlread('d.stl');

V = model.Points; % 顶点坐标

F = model.ConnectivityList; % 面索引（三角面）

Z = V(:,3); % 每个顶点的 Z 高度作为颜色依据

figure;

patch('Faces', F, ...
'Vertices', V, ...
'FaceVertexCData', Z, ...
'FaceColor', 'interp', ... % 插值显示
'EdgeColor', 'k', ...
'FaceAlpha', 0.5, ...
'LineWidth', 0.3);
colormap(turbo); % 渐变色，可换为 parula、jet 等
% colorbar;
axis equal;
view(3);
axis off;          % 隐藏坐标轴、刻度、边框
camzoom(1.9);  % 放大2倍，可以调整这个数值，比如 1.5, 2, 3...
% xlabel('X'); ylabel('Y'); zlabel('Z');

% title('按高度渐变显示 STL 模型');

% ===== 显示所有顶点 =====

hold on; % 重要：确保后续绘图不会覆盖 patch

% plot3(Vertices(:,1), Vertices(:,2), Vertices(:,3), ...
% 'bo', 'MarkerSize', 1, 'MarkerFaceColor', 'b'); % 蓝色实心圆点

% ===== 显示所有顶点的法向量（红色箭头）=====

% 为了避免太密集，可以选择显示全部，或每隔N个点显示一个

N_display = 1; % 每隔N个顶点显示一个法向量，可调

sample_verts_idx = 1:N_display:size(Vertices,1); % 采样顶点索引

verts_to_show = Vertices(sample_verts_idx, :); % 采样顶点位置

normals_to_show = vertex_normals(sample_verts_idx, :); % 对应法向量

% 设置箭头长度（可调，原始法向量可能太短）

normal_length_scale = 0.01;

quiver3(verts_to_show(:,1), verts_to_show(:,2), verts_to_show(:,3), ...
normals_to_show(:,1)*normal_length_scale, ...
normals_to_show(:,2)*normal_length_scale, ...
normals_to_show(:,3)*normal_length_scale, ...
'b', 'LineWidth', 1.2, 'MaxHeadSize', 0.5, 'AutoScale', 'off');
axis equal;


hold off;