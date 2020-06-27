pMid = 0.5 * (pMesh(:, N / 2) + pMesh(:, N / 2 + 1));
uMagMid = 0.5 * (UMagMesh(:, N / 2) + UMagMesh(:, N / 2 + 1));
uMid = 0.5 * (UMesh(:, N / 2) + UMesh(:, N / 2 + 1));
yMid = y;

yMid = [0, y, 1];
uMid = [0; uMid; 1];
yMidTest = [0, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531, 0.5, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 0.9766, 1];
uMidTest = [0, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662, -0.21090, -0.20581, -0.13641, 0.00332, 0.23151, 0.68717, 0.73772, 0.78871, 0.84123, 1];

% y_fluent = y_fluent + 0.5;

% figure(4)
% hold on
% plot(yMid, pMid);
% % plot(y_fluent + 0.5, p_fluent);
% hold off

figure(5)
hold on
plot(yMid, uMid, 'r-', 'Linewidth', 2);
plot(yMidTest, uMidTest, 'sk', 'markerFaceColor', 'k')
% plot(y_fluent + 0.5, magU_fluent);
xlabel('y / m');
ylabel('u / (m/s)');
legend('My results', 'Benchmark results');
hold off