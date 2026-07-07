function M = dvh_metrics(dose, Rx, Vlevels)
% DVH_METRICS  DVH metrics for a single structure.
%
%   M = dvh_metrics(dose)              % no prescription: V95/HI not computed
%   M = dvh_metrics(dose, Rx)          % with prescription Rx (Gy)
%   M = dvh_metrics(dose, Rx, Vlevels) % + Vx at the dose levels in Vlevels (Gy)
%
%   dose    : per-voxel dose (Gy) of the structure.
%   Rx      : reference prescription dose (Gy). Use for the target.
%   Vlevels : dose levels (Gy) for Vx = % of the volume receiving >= x Gy.
%             Ex.: [15 20 25]. Accessible via M.V: M.V.V15, M.V.V20, M.V.V25.
%
%   ICRU-83 convention:  Dx = dose received by x% of the volume.
%   Fields: Dmean, Dmax, Dmin, D95, D98, D2, D50, V95, HI [, V.Vxx ...]
%
%   Toolbox-free (percentile implemented locally).
    dose = dose(:);
    M.Dmean = mean(dose);
    M.Dmax  = max(dose);
    M.Dmin  = min(dose);
    M.D95 = pctl(dose,  5);   % dose covering 95% of the volume
    M.D98 = pctl(dose,  2);   % near-minimum (ICRU-83)
    M.D2  = pctl(dose, 98);   % near-maximum (ICRU-83)
    M.D50 = pctl(dose, 50);   % median
    if nargin > 1 && ~isempty(Rx)
        M.V95 = 100 * mean(dose >= 0.95*Rx);   % % of volume with >= 95% of Rx
        M.HI  = (M.D2 - M.D98) / M.D50;        % homogeneity index (ICRU-83)
    else
        M.V95 = NaN;
        M.HI  = NaN;
    end
    % generic Vx (percent of volume receiving >= level Gy)
    if nargin > 2 && ~isempty(Vlevels)
        for lv = Vlevels(:)'
            campo = sprintf('V%g', lv);
            M.V.(campo) = 100 * mean(dose >= lv);
        end
    end
end
function q = pctl(v, p)
% Percentile by linear interpolation (same method as MATLAB's prctile),
% without requiring the Statistics Toolbox.
    v = sort(v(:));
    n = numel(v);
    if n == 1, q = v; return; end
    idx  = (p/100)*(n-1) + 1;
    lo   = floor(idx);
    hi   = ceil(idx);
    frac = idx - lo;
    q = v(lo)*(1-frac) + v(hi)*frac;
end