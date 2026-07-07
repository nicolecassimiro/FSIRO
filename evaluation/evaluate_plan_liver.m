function T = evaluate_plan_liver(dose, sizes, Dp, names)
% EVALUATE_PLAN_LIVER  Clinical (DVH) evaluation of a liver SBRT plan.
%
%   T = evaluate_plan_liver(out.dose, out.sizes, out.Dp, out.names)
%
% Robust to variable structure sets: each structure present in the plan is
% matched to its QUANTEC spec by name, so patients missing an OAR (e.g. no
% pancreas) or with extra ones (e.g. bowels) are handled automatically. If
% NAMES is omitted, the standard protocol_liver order is assumed.
%
% Uses DVH_METRICS as the engine. QUANTEC limits are conventional-fractionation
% reference values (SBRT limits differ; audit for the fractionation scheme).
% Returns T (struct array), one entry per structure.

    if nargin < 4 || isempty(names)
        P = protocol_liver();
        names = [{P.target.pattern}, {P.oars.pattern}];
    end
    names = cellfun(@clean_name, names, 'uni', 0);

    % QUANTEC spec map: key = structure name -> {metric, limit, type, label}
    spec = containers.Map('KeyType','char','ValueType','any');
    spec('Liver minus CTV') = {'Dmean', 28,  'parallel', 'RILD (primary/HCC)'};
    spec('Liver minus GTV') = {'Dmean', 28,  'parallel', 'RILD (primary/HCC)'};
    spec('Spinal Cord')     = {'Dmax',  50,  'serial',   'myelopathy (SBRT ref ~18)'};
    spec('Heart')           = {'Dmean', 26,  'parallel', 'pericarditis (aim V25<10%)'};
    spec('Stomach')         = {'Dmin',  45,  'serial',   'D100 whole-stomach'};
    spec('Oesophagus')      = {'Dmean', 34,  'parallel', 'esophagitis'};
    spec('Duodenum')        = {'Dmax',  45,  'serial',   'small bowel'};
    spec('Bowels')          = {'Dmax',  45,  'serial',   'small bowel (V15<120cc)'};
    spec('Pancreas')        = {'Dmean', NaN, 'parallel', 'no specific QUANTEC'};
    spec('Kidney (L)')      = {'Dmean', 18,  'parallel', 'renal dysfunction'};
    spec('Kidney (R)')      = {'Dmean', 18,  'parallel', 'renal dysfunction'};

    % slice the dose per structure
    D = cell(1,numel(sizes));  o = 0;
    for i = 1:numel(sizes), D{i} = dose(o+1:o+sizes(i));  o = o + sizes(i); end

    fprintf('=== Liver (SBRT) clinical evaluation  (Dp = %.1f Gy) ===\n', Dp);

    % ---------------- target (first structure) ----------------
    Mp = dvh_metrics(D{1}, Dp);
    covmin = 0.67 * Mp.Dmax;                 % SBRT: prescribed minimum = 67% of max
    okD98  = Mp.D98 >= covmin;
    fprintf(['%-16s: D98=%.2f (>=%.2f? %s)  D2=%.2f  D50=%.2f  Dmean=%.2f  ' ...
             'Dmax=%.2f  HI=%.3f\n'], names{1}, Mp.D98, covmin, yn(okD98), ...
             Mp.D2, Mp.D50, Mp.Dmean, Mp.Dmax, Mp.HI);
    T(1) = res(names{1}, 'D98', Mp.D98, covmin, okD98);

    % ---------------- OARs (match each present structure by name) ----------------
    for i = 2:numel(sizes)
        nm = names{i};  M = dvh_metrics(D{i});
        if isKey(spec, nm)
            s = spec(nm);  fld = s{1};  lim = s{2};  typ = s{3};  lab = s{4};
        else
            fld = 'Dmean';  lim = NaN;  typ = 'parallel';  lab = 'no spec';
        end
        val = M.(fld);
        if isnan(lim)
            fprintf('%-16s: %s=%.2f Gy  [%s]  (%s)\n', nm, fld, val, typ, lab);
            ok = true;
        else
            ok = val <= lim;
            fprintf('%-16s: %s=%.2f Gy (<=%g? %s)  [%s]  (%s)\n', ...
                    nm, fld, val, lim, yn(ok), typ, lab);
        end
        T(end+1) = res(nm, fld, val, lim, ok); %#ok<AGROW>
    end
end

function s = yn(tf)
    if tf, s = 'yes'; else, s = 'NO'; end
end

function r = res(name, metric, value, limit, passed)
    r.name = name;  r.metric = metric;  r.value = value;
    r.limit = limit;  r.passed = passed;
end

function s = clean_name(pat)
    s = regexprep(pat, '[\^\$]', '');
    s = regexprep(s, '\((CTV\|GTV)\)', 'CTV');   % '(CTV|GTV)' -> 'CTV' for display
    s = strrep(strrep(s, '\(', '('), '\)', ')');
end