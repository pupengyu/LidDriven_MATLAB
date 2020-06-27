%% 后处理
figure(1)
hold on
contourf(X, Y, pMesh, 1000, 'edgeColor', 'none');
colormap(jet(1000));
caxis([-0.2, 0.4]);
colorbar;

% q = quiver(X, Y, UMesh, VMesh);
% q.AutoScale = 'on';
% axis([0, 1, 0, 1]);
% axis square;
hold off

figure(2)
q = quiver(X, Y, UMesh, VMesh);
q.AutoScale = 'on';
axis([0, 1, 0, 1]);
axis square;

figure(3)
hold on;
contourf(X, Y, UMagMesh, 1000, 'edgeColor', 'none');
colormap(jet(1000));
caxis([0, 1]);
colorbar;

q = quiver(X, Y, UMesh, VMesh, 5, 'w', 'fill');
q.AutoScale = 'on';
axis([0, 1, 0, 1]);
axis square;
hold off