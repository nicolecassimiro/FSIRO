function P = protocol_liver()
% PROTOCOL_LIVER  Liver SBRT protocol for FSIRO (full OAR set).
% dT (=75) and Dp are read from the .mat automatically.
%
% FEASIBILITY: the target is stereotactic (heterogeneous, minimum ~67% of the
% maximum), so the target membership keeps a high peak but a feasible lower
% support. Each OAR 'support' (bE) is set ABOVE the achievable per-voxel MAX
% dose (measured on a converged plan), otherwise the omega formulation becomes
% infeasible (Ax <= bE cannot hold for voxels adjacent to the PTV). In
% particular the healthy liver touches the PTV and receives up to ~dT, so its
% support is ~dT. Clinical dose limits are applied in the DVH evaluation, not
% as supports; the low 'goal' drives sparing.

    P.name   = 'Liver';
    P.mu_min = 0.2;                % no P.constraints -> pure fuzzy formulation

    % Target: high peak toward the 75 Gy maximum, feasible lower support.
    P.target.pattern = '^PTV$';
    P.target.exclude = '';
    P.target.lo      = 0.60;
    P.target.core_lo = 0.90;
    P.target.core_hi = 0.98;
    P.target.peak    = 1.000;
    P.target.hi      = 1.070;

    % OARs: support > measured per-voxel MAX; goal (plateau) low for sparing.
    % (measured MAX on converged plan, Gy): liver 74.5, kidneyR 49.7, oeso 30.7,
    %  heart 25.8, pancreas 21.2, stomach 16.4, duodenum 14.4, kidneyL 10.2,
    %  cord 9.1 -> supports set safely above each.
    o = @(pat,g,sup,ty) struct('pattern',pat,'goal',g,'support',sup,'type',ty);
    P.oars = [ ...
        o('^Liver minus CTV',10,  78, 'parallel'), ...
        o('^Spinal Cord',    10,  40, 'serial'  ), ...
        o('^Heart',          10,  50, 'serial'  ), ...
        o('^Stomach',        10,  40, 'serial'  ), ...
        o('^Oesophagus',     10,  45, 'serial'  ), ...
        o('^Duodenum',       10,  40, 'serial'  ), ...
        o('^Pancreas',       10,  40, 'serial'  ), ...
        o('^Kidney \(L\)',   10,  40, 'parallel'), ...
        o('^Kidney \(R\)',   10,  55, 'parallel') ];
end