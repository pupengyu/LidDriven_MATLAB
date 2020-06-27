fp = fopen('result.dat', 'w');

fprintf(fp, 'title = \"result\"\n');
fprintf(fp, 'variables = x, y, p, u, v, magU\n');
fprintf(fp, 'ZONE I = 10, J = 10, F = Point\n');
x = reshape(X', N * N, 1);
y = reshape(Y', N * N, 1);
magU = reshape(UMagMesh', N * N, 1);
results = [x, y, pNew, UNew(:, 1), UNew(:, 2), magU];
for ii = 1:N * N
    for jj = 1:6
        fprintf(fp, '%f\t', results(ii, jj));
    end
    fprintf(fp, '\n');
end

fclose(fp);