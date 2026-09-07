% =============================================
% MATLAB 脚本：从 STL 模型表面均匀采样点
% 方法：基于三角面片面积加权随机采样
% =============================================

% 1. 读取 STL 文件
filename = 'd.stl'; % 替换为你的 STL 文件路径
stl = stlread(filename);     % 需要 stlread 函数（见下文说明）

% 2. 提取顶点和面片
vertices = stl.Points;       % Nx3 矩阵，所有顶点 [x,y,z]
faces = stl.ConnectivityList;% Mx3 矩阵，每个三角形由三个顶点索引组成

% 3. 计算每个三角形的面积
numFaces = size(faces, 1);
areas = zeros(numFaces, 1);

for i = 1:numFaces
    % 获取当前三角形的三个顶点
    v1 = vertices(faces(i,1), :);
    v2 = vertices(faces(i,2), :);
    v3 = vertices(faces(i,3), :);
    
    % 计算两个边向量
    vec1 = v2 - v1;
    vec2 = v3 - v1;
    
    % 叉积求面积（三角形面积为叉积模长的一半）
    cross_prod = cross(vec1, vec2);
    area = 0.5 * norm(cross_prod);
    areas(i) = area;
end

% 4. 构造按面积加权的概率分布
totalArea = sum(areas);
probabilities = areas / totalArea;  % 每个面片被选中的概率

% 5. 设定要采样的点的总数
numSamples = 15000; % 你可以修改为你需要的点数，比如 5000

% 6. 按概率随机选择面片索引
selectedFaceIndices = randsample(1:numFaces, numSamples, true, probabilities);

% 7. 在每个选中的三角面上随机均匀采样一个点
sampledPoints = zeros(numSamples, 3);

for i = 1:numSamples
    faceIdx = selectedFaceIndices(i);
    v1 = vertices(faces(faceIdx, 1), :);
    v2 = vertices(faces(faceIdx, 2), :);
    v3 = vertices(faces(faceIdx, 3), :);
    
    % 在三角面片上均匀随机采样一个点
    % 方法：重心坐标法
    r1 = sqrt(rand());  % 保证在三角形内均匀分布，不是单纯随机
    r2 = rand();
    
    if (r1 + r2 > 1)
        r1 = 1 - r1;
        r2 = 1 - r2;
    end
    
    % 重心坐标公式
    p = (1 - r1 - r2) * v1 + r1 * v2 + r2 * v3;
    sampledPoints(i, :) = p;
end

% 8. （可选）可视化结果
figure;
hold on;
trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3), ...
    'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); % 原始模型
scatter3(sampledPoints(:,1), sampledPoints(:,2), sampledPoints(:,3), ...
    10, 'r', 'filled'); % 采样点，红色
axis equal;
title(['在 STL 表面均匀采样 ' num2str(numSamples) ' 个点']);
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;
hold off;

% 9. （可选）导出采样点
% sampledPoints 是一个 N×3 矩阵，可以直接用于其它用途
% 例如保存为 .mat 或 .csv
% save('sampled_points.mat', 'sampledPoints');
% writematrix(sampledPoints, 'sampled_points.csv');