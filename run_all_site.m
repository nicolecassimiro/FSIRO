function run_all_site(cfg)
% RUN_ALL_SITE  Batch-run FSIRO over all patients of one site, both membership
% shapes, and record DVH metrics for triangular, trapezoidal and TROTS.
%
% Saves, per plan, a .mat under <outdir>/planos/ (x, dose, metrics, ...), and
% appends rows to <outdir>/<site>_resumo.csv (long/tidy: one row per
% patient x structure x method). Resumable: plans whose .mat already exists are
% not re-optimized. TROTS enters as the method 'TROTS' (reference column).
%
% cfg fields: site, datadir, prefix, npat, protocol (handle), Vlevels, outdir.

    P = cfg.protocol();
    labels = [{clean_name(P.target.pattern)}, cellfun(@clean_name, {P.oars.pattern}, 'uni', 0)];

    planos_dir = fullfile(cfg.outdir, 'planos');
    if ~exist(planos_dir, 'dir'), mkdir(planos_dir); end
    csvfile = fullfile(cfg.outdir, [cfg.site '_resumo.csv']);

    Vcols = arrayfun(@(v) sprintf('V%g', v), cfg.Vlevels, 'uni', 0);
    cols  = [{'patient','structure','method','D98','D95','D2','D50','Dmean', ...
              'Dmax','Dmin','V95','HI'}, Vcols];
    if ~exist(csvfile, 'file')
        fid = fopen(csvfile, 'w');  fprintf(fid, '%s\n', strjoin(cols, ','));  fclose(fid);
    end

    for k = 1:cfg.npat
        pid  = sprintf('%s_%02d', cfg.prefix, k);
        mfile = fullfile(cfg.datadir, [pid '.mat']);
        if ~exist(mfile, 'file'), warning('missing %s -- skipping', mfile); continue; end
        fprintf('==== %s (%d/%d) ====\n', pid, k, cfg.npat);

        outs = struct();  ok = true;
        for f = {'triangular','trapezoidal'}
            shape = f{1};
            planmat = fullfile(planos_dir, sprintf('%s_%02d_%s.mat', cfg.site, k, shape));
            if exist(planmat, 'file')
                S = load(planmat);  out = S.out;
                fprintf('   %-12s: cached\n', shape);
            else
                try
                    t = tic;  out = FSIRO(mfile, shape, P, false);
                    fprintf('   %-12s: %d it, %.0f s, cp=%.1e (%s)\n', ...
                            shape, out.iter, out.time, out.cp, out.criterion);
                    save(planmat, 'out');
                catch ME
                    warning('FSIRO failed on %s %s: %s', pid, shape, ME.message);
                    ok = false;  break;
                end
            end
            outs.(shape) = out;
        end
        if ~ok, continue; end

        sizes = outs.triangular.sizes;  Dp = outs.triangular.Dp;
        dref  = trots_dose(mfile, P);   % TROTS reference, same structure order

        fid = fopen(csvfile, 'a');
        rid = sprintf('%s_%02d', cfg.site, k);
        append_rows(fid, rid, 'triangular',  outs.triangular.dose,  sizes, labels, Dp, cfg.Vlevels);
        append_rows(fid, rid, 'trapezoidal', outs.trapezoidal.dose, sizes, labels, Dp, cfg.Vlevels);
        if ~isempty(dref)
            append_rows(fid, rid, 'TROTS',    dref,                 sizes, labels, Dp, cfg.Vlevels);
        end
        fclose(fid);
    end
    fprintf('\nDone. Plans: %s   Summary: %s\n', planos_dir, csvfile);
end


% ---------- append one row per structure to the CSV ----------
function append_rows(fid, patient, method, dose, sizes, labels, Dp, Vlevels)
    o = 0;
    for i = 1:numel(sizes)
        d = dose(o+1 : o+sizes(i));  o = o + sizes(i);
        if i == 1, M = dvh_metrics(d, Dp, Vlevels);   % PTV: with prescription
        else,      M = dvh_metrics(d, [], Vlevels);   % OAR: no prescription
        end
        fprintf(fid, '%s,%s,%s,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g', ...
                patient, labels{i}, method, M.D98, M.D95, M.D2, M.D50, ...
                M.Dmean, M.Dmax, M.Dmin, M.V95, M.HI);
        for lv = Vlevels(:)'
            fprintf(fid, ',%.4g', M.V.(sprintf('V%g', lv)));
        end
        fprintf(fid, '\n');
    end
end


% ---------- TROTS reference dose, structures in protocol order ----------
function dref = trots_dose(mfile, P)
    S = load(mfile);
    if ~isfield(S, 'solutionX') || isempty(S.solutionX), dref = []; return; end
    data = S.data;
    pats = [{P.target.pattern}, {P.oars.pattern}];
    A = [];
    for i = 1:numel(pats)
        idx = find(~cellfun('isempty', regexpi({data.matrix.Name}, pats{i}, 'once')) & ...
                    cellfun('isempty', regexpi({data.matrix.Name}, '\(mean\)', 'once')), 1);
        A = [A; double(data.matrix(idx).A)];
    end
    dref = A * S.solutionX;
end


% ---------- turn a regex pattern into a readable structure label ----------
function s = clean_name(pat)
    s = regexprep(pat, '[\^\$]', '');
    s = strrep(strrep(s, '\(', '('), '\)', ')');
end