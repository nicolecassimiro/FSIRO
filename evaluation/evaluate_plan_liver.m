function T = evaluate_plan_liver(dose, sizes, Dp)
% EVALUATE_PLAN_LIVER  Clinical (DVH) evaluation of a liver SBRT plan.
%
%   T = evaluate_plan_liver(out.dose, out.sizes, out.Dp)
%
%   Uses METRICAS_DVH as the engine. Structure order follows protocol_liver:
%     PTV, Liver minus CTV, Spinal Cord, Heart, Stomach, Oesophagus, Duodenum,
%     Pancreas, Kidney (L), Kidney (R).
%
%   QUANTEC limits (conventional fractionation; see criterios_dose_Liver). SBRT
%   limits differ - the values below are reference and should be audited for the
%   fractionation scheme. Returns T (struct array), one entry per structure.

    names = {'PTV','Liver minus CTV','Spinal Cord','Heart','Stomach', ...
             'Oesophagus','Duodenum','Pancreas','Kidney (L)','Kidney (R)'};

    % slice the dose per structure
    D = cell(1,numel(sizes));  o = 0;
    for i = 1:numel(sizes), D{i} = dose(o+1:o+sizes(i));  o = o + sizes(i); end

    fprintf('=== Liver (SBRT) clinical evaluation  (Dp = %.1f Gy) ===\n', Dp);

    % ---------------- PTV (target) ----------------
    Mp = dvh_metrics(D{1}, Dp);
    covmin = 0.67 * Mp.Dmax;                 % SBRT: prescribed minimum = 67% of max
    okD98  = Mp.D98 >= covmin;
    fprintf(['PTV   : D98=%.2f (>=%.2f? %s)  D2=%.2f  D50=%.2f  Dmean=%.2f  ' ...
             'Dmax=%.2f  HI=%.3f\n'], Mp.D98, covmin, yn(okD98), Mp.D2, Mp.D50, ...
             Mp.Dmean, Mp.Dmax, Mp.HI);
    T(1) = res('PTV','D98',Mp.D98,covmin,okD98);

    % ---------------- OARs ----------------
    % spec: {index, metric, limit, type, label}
    spec = { ...
        2, 'Dmean', 28,  'parallel', 'RILD (primary/HCC)'; ...
        3, 'Dmax',  50,  'serial',   'myelopathy (SBRT ref ~18)'; ...
        4, 'Dmean', 26,  'parallel', 'pericarditis (aim V25<10%)'; ...
        5, 'Dmin',  45,  'serial',   'D100 whole-stomach'; ...
        6, 'Dmean', 34,  'parallel', 'esophagitis'; ...
        7, 'Dmax',  45,  'serial',   'small bowel'; ...
        8, 'Dmean', NaN, 'parallel', 'no specific QUANTEC'; ...
        9, 'Dmean', 18,  'parallel', 'renal dysfunction'; ...
       10, 'Dmean', 18,  'parallel', 'renal dysfunction' };

    for k = 1:size(spec,1)
        idx = spec{k,1};  fld = spec{k,2};  lim = spec{k,3};
        typ = spec{k,4};  lab = spec{k,5};
        M   = dvh_metrics(D{idx});
        val = M.(fld);
        if isnan(lim)
            fprintf('%-15s: %s=%.2f Gy  [%s]  (%s)\n', names{idx}, fld, val, typ, lab);
            ok = true;
        else
            ok = val <= lim;
            fprintf('%-15s: %s=%.2f Gy (<=%d? %s)  [%s]  (%s)\n', ...
                    names{idx}, fld, val, lim, yn(ok), typ, lab);
        end
        T(end+1) = res(names{idx}, fld, val, lim, ok); %#ok<AGROW>
    end
end

function s = yn(tf)
    if tf, s = 'yes'; else, s = 'NO'; end
end

function r = res(name, metric, value, limit, passed)
    r.name = name;  r.metric = metric;  r.value = value;
    r.limit = limit;  r.passed = passed;
end