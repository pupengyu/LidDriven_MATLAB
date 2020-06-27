%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 这是一个用于计算Lid driven flow的Matlab程序，采用的是有限体积法，%
% SIMPLE算法，对流项一阶迎风格式，扩散项中心差分格式，界面质量流量采 %
% 用Rhie-Chow插值。使用Collocated二维结构化网格。计算工况的雷诺数为 %
% 100                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a Matlab code for the lid driven flow simulation using %
% finite volume method with SIMPLE algorithm. The convection terms %
% are discretized by the first order upwind scheme, and the diffusion  %
% term the central difference scheme. A set of collocated structured %
% 2D mesh is used with the Rhie-Chow interpolation. Re = 100 %                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 清除数据与图像
clear
close all
clc

%% 初始化
tic     % 计时起点
profile on;
global L mu h ub alphaU alphaP
L = 1.0;    % Lid长度
mu = 0.01;  % 动力粘度系数
ub = 1;     % Lid运动速度
alphaU = 0.7;   % 动量方程隐式松弛因子
alphaP = 0.3;   % 压力修正显式松弛因子
Error = 1e-5;   % 收敛准则

N = 32; % 一条边上的网格数，总网格数为N^2
h = L / N;  % 网格尺寸

UOld = zeros(N * N, 2); % 上一个迭代步中的速度场
UNew = zeros(N * N, 2); % 该迭代步中的速度场
p = zeros(N * N, 1);    % 上一个迭代步中的压强场
pNew = zeros(N * N, 1); % 该迭代步中的压强场
pC = zeros(N * N, 1);   % 压力修正

v_RES = zeros(100000, 1);   % 残差序列，用于输出残差图

%% 设置单元类型，用以识别边界条件
cellType = cell(N * N, 1);
for ii = 1:N*N
    cellType{ii} = getCellType(ii, N);
end

%% SIMPLE求解
for ii = 1:100000
    %% 求解动量方程，返回满足动量方程的速度场，与中间变量D，后者由于Rhie-Chow插值
    [UNew, D] = predictU(cellType, UOld, p, N);
    
    %% 求解压力修正方程，并修正压强与速度场
    [pNew, UNew] = correctP(cellType, UNew, p, D, N);

    %% 获取残差并判断是否收敛
    RES = max(max(abs(UOld - UNew)));   % 这里用的是error而非residual，主要是嫌麻烦
    v_RES(ii) = RES;
    if RES < Error
        disp('Converged!');
        break; 
    else
        disp(['RES = ', num2str(RES)]);
        figure(1);
        plot(1:ii, log10(v_RES(1:ii)), '-');    % 用于显示残差数值并绘制残差曲线，会影响计算速度，可以注释掉
    end
    
    %% 更新物理场
    UOld = UNew;
    p = pNew;
end

toc     % 计时终点
myProfile = profile('info');

%% 后处理
x = h * (1:N) - 0.5 * h;
y = x;
[X, Y] = meshgrid(x, y);
UMesh = reshape(UOld(:, 1), N, N)';
VMesh = reshape(UOld(:, 2), N, N)';
UMagMesh = sqrt(UMesh .^ 2 + VMesh .^ 2);
pMesh = reshape(p, N, N)';
[sx, sy] = meshgrid(0:0.25:1, 0:0.1:1);

showPlot;
exportMidLineValue;

% % 压强云图
% figure(1)
% contourf(X, Y, pMesh, 100, 'edgeColor', 'none');
% colormap(jet(100));
% caxis([-0.2, 0.4]);
% colorbar;
% 
% % 速度矢量图
% figure(2)
% q = quiver(X, Y, UMesh, VMesh);
% q.AutoScale = 'on';
% axis([0, 1, 0, 1]);
% axis square;
% 
% % 速度大小云图
% figure(3)
% contourf(X, Y, UMagMesh, 100, 'edgeColor', 'none');
% colormap(jet(100));
% caxis([0, 1]);
% colorbar;