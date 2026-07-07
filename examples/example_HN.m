% EXAMPLE_HN  Minimal usage of FSIRO on a TROTS Head-and-Neck case.
% Adjust MATFILE to point at your TROTS .mat file.

matfile = fullfile('Head-and-Neck', 'Head-and-Neck_01.mat');
P = protocol_HN();

out_tri  = FSIRO(matfile, 'triangular',  P);
out_trap = FSIRO(matfile, 'trapezoidal', P);

fprintf('\n%-14s | %10s | %10s\n', 'Structure', 'Triangular', 'Trapezoidal');
labels = [{P.target.pattern}, {P.oars.pattern}];
for i = 1:numel(out_tri.means)
    fprintf('%-14s | %10.2f | %10.2f\n', labels{i}, out_tri.means(i), out_trap.means(i));
end
fprintf('\nStops: %s (%s) | %s (%s)\n', ...
        out_tri.criterion, out_tri.shape, out_trap.criterion, out_trap.shape);