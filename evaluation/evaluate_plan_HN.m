function T = evaluate_plan_HN(dose, sizes, Dp)
% EVALUATE_PLAN_HN  Clinical (DVH) evaluation of a head-and-neck plan for all
% structures in protocol_HN, applying ICRU-83 to the target and QUANTEC limits
% to the OARs (serial -> Dmax, parallel -> Dmean).
%
%   T = evaluate_plan_HN(out.dose, out.sizes, out.Dp)
%
% Uses DVH_METRICS as the engine and reads structure names, types and clinical
% limits from protocol_HN. Returns T (struct array), one entry per structure.

    P = protocol_HN();
    labels = [{clean_name(P.target.pattern)}, cellfun(@clean_name, {P.oars.pattern}, 'uni', 0)];

    % slice the dose per structure (target first, then OARs)
    D = cell(1,numel(sizes));  o = 0;
    for i = 1:numel(sizes), D{i} = dose(o+1:o+sizes(i));  o = o + sizes(i); end

    fprintf('=== Head-and-Neck clinical evaluation  (Dp = %.1f Gy) ===\n', Dp);

    % ---------------- PTV (ICRU-83) ----------------
    Mp = dvh_metrics(D{1}, Dp);
    okD98 = Mp.D98 >= 0.95*Dp;   okD2 = Mp.D2 <= 1.07*Dp;
    fprintf(['PTV            : D98=%.2f (>=%.2f? %s)  D2=%.2f (<=%.2f? %s)  ' ...
             'D50=%.2f  Dmean=%.2f  V95=%.1f%%  HI=%.3f\n'], ...
             Mp.D98, 0.95*Dp, yn(okD98), Mp.D2, 1.07*Dp, yn(okD2), ...
             Mp.D50, Mp.Dmean, Mp.V95, Mp.HI);
    T(1) = res(labels{1}, 'D98', Mp.D98, 0.95*Dp, okD98);

    % ---------------- OARs ----------------
    for j = 1:numel(P.oars)
        oar = P.oars(j);  M = dvh_metrics(D{j+1});
        if strcmp(oar.type,'serial'), fld = 'Dmax'; else, fld = 'Dmean'; end
        val = M.(fld);  lim = oar.eval_limit;  ok = val <= lim;
        fprintf('%-15s: %s=%.2f Gy (<=%g? %s)  [%s]\n', ...
                labels{j+1}, fld, val, lim, yn(ok), oar.type);
        T(end+1) = res(labels{j+1}, fld, val, lim, ok); %#ok<AGROW>
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
    s = strrep(strrep(s, '\(', '('), '\)', ')');
end