function out = FSIRO(matfile, shape, protocol, verbose)
% FSIRO  Fuzzy Surprise Interior-point Radiotherapy Optimization.
% Beamlet-intensity optimization of a fuzzy-surprise dose model, solved by a
% primal-dual interior-point method with Mehrotra's predictor-corrector.
%
%   out = FSIRO(matfile, shape)                    % Head-and-Neck (default)
%   out = FSIRO(matfile, shape, protocol)          % custom protocol
%   out = FSIRO(matfile, shape, protocol, true)    % print each iteration
%
% Model (fuzzy support enforced as a hard constraint):
%     min_x  sum_i S(mu_i(A x))   s.t.  A x = w,  b1 <= w <= bE,  0 <= x <= U
% where S is Neumaier's surprise (1/mu - 1)^2, regularized to be C^1 (see
% MEMBERSHIP_SMOOTH), and mu is a fuzzy membership whose shape is 'triangular'
% (3-point) or 'trapezoidal' (4-point). Structures are located by name from the
% protocol, so no per-patient index editing is required.
%
% INPUT
%   matfile  : path to a TROTS-format .mat file (variables data, problem).
%   shape    : 'triangular' | 'trapezoidal'.
%   protocol : struct from protocol_HN() or an equivalent (see that file).
%   verbose  : true to print per-iteration progress (default false).
%
% OUTPUT (struct)
%   x, dose, phi, iter, time, cp, criterion, means, sizes, names.
%
% The Newton step is solved matrix-free by preconditioned CG; the method stops
% when the plan stabilizes (mean structure doses change < tol_plan over a
% window) or the KKT residual falls below tol_kkt.

    if nargin < 3 || isempty(protocol), protocol = protocol_HN(); end
    if nargin < 4 || isempty(verbose),  verbose  = false; end
    shape = lower(shape);
    tri = strcmp(shape,'triangular');
    if ~tri && ~strcmp(shape,'trapezoidal')
        error('shape must be ''triangular'' or ''trapezoidal''.');
    end

    tau      = 0.99995;   % fraction-to-the-boundary
    tol_kkt  = 1e-4;      % KKT stopping tolerance
    tol_cg   = 1e-6;      % inner CG tolerance
    tol_plan = 0.05;      % plan-stability tolerance (Gy)
    win      = 5;         % consecutive stable iterations to stop
    maxit    = 300;
    mu_min   = protocol.mu_min;

    % ---- load and locate structures by name (target first, then OARs) ----
    Sm = load(matfile);  data = Sm.data;  problem = Sm.problem;
    At = full_matrix(data, protocol.target.pattern, protocol.target.exclude);
    nOAR = numel(protocol.oars);
    Aoar = cell(1,nOAR);  szoar = zeros(1,nOAR);  names = cell(1,nOAR+1);
    names{1} = protocol.target.pattern;
    for j = 1:nOAR
        Aoar{j} = full_matrix(data, protocol.oars(j).pattern, '');
        szoar(j) = size(Aoar{j},1);
        names{j+1} = protocol.oars(j).pattern;
    end
    A = [At; vertcat(Aoar{:})];  A = sparse(A);  A2 = A.^2;
    [TP, N] = size(A);
    szt = size(At,1);  sizes = [szt szoar];
    U = 1e3*ones(N,1);

    % ---- target reference dose dT and prescription Dp from the database ----
    [dT, Dp] = target_scale(problem, protocol.target, protocol);

    % ---- fuzzy numbers b (per structure), shape-dependent ----
    b = build_b(protocol, shape, szt, szoar, dT);
    bL = b(:,1);  bE = b(:,end);

    % ---- strictly interior starting point ----
    x = init_point(data);  x = min(max(x,1), U-1);
    w = A*x;  w = min(max(w, bL+1e-3*(bE-bL)), bE-1e-3*(bE-bL));
    p = 1e-2*ones(N,1);  v = U - x;  s = 10*ones(N,1);
    z1 = w - bL;  z2 = bE - w;  t1 = 10*ones(TP,1);  t2 = 10*ones(TP,1);
    y = 200*ones(N,1);  w1 = 20*ones(TP,1);  w2 = -20*ones(TP,1);  q = A*x - w;

    % ---- optional linear inequality constraints  G x >= 0 (PTV min/max dose) ----
    use_con = isfield(protocol,'constraints') && protocol.constraints;
    if use_con, G = find_constraints(data, problem, N); else, G = sparse(0,N); end
    mG = size(G,1);  G2 = G.^2;
    h   = max(G*x, 1);            % slack h = G x > 0 (strictly interior)
    eta = ones(mG,1);            % multiplier of G x >= 0

    npair = 2*N + 2*TP + mG;
    means_prev = inf(1,nOAR+1);  stable = 0;  cp = 1;  it = 0;
    criterion = sprintf('maxit=%d', maxit);
    t0 = clock;

    while it < maxit
        d = A*x;
        [sup, dg, dH] = membership_smooth(d, b, mu_min);
        gradf = A'*dg;

        r1 = A'*q + y + p + G'*eta - gradf;
        r2 = w1 + w2 - q;
        r3 = -s - y;
        r4 = t1 - w1;
        r5 = -t2 - w2;
        r6 = U - x - v;
        r7 = bL + z1 - w;
        r8 = bE - z2 - w;
        r9 = w - A*x;
        rH = G*x - h;                       % slack def:  G x - h = 0

        xp = x.*p;  vs = v.*s;  z1t1 = z1.*t1;  z2t2 = z2.*t2;  he = h.*eta;
        Dd = p./x + s./v;  D = dH + t1./z1 + t2./z2;  ge = eta./h;

        % predictor (affine, no centering)
        r10=-xp; r11=-vs; r12=-z1t1; r13=-z2t2; rHE=-he;
        rc = reduce_rhs(A,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,s,v,t1,t2,z1,z2,x) ...
             + G'*((rHE + eta.*rH)./h);
        dxa = cg_step(A, A2, D, Dd, G, G2, ge, rc, tol_cg);
        wa  = A*dxa - r9;  dz1a = wa - r7;  dz2a = r8 - wa;
        dva = r6 - dxa;    dsa  = (r11 - s.*dva)./v;
        dpa = (r10 - p.*dxa)./x;
        dt1a = (r12 - t1.*dz1a)./z1;  dt2a = (r13 - t2.*dz2a)./z2;
        dha  = G*dxa - rH;   detaa = (rHE - eta.*dha)./h;
        apa = min([bstep(x,dxa,1) bstep(v,dva,1) bstep(z1,dz1a,1) bstep(z2,dz2a,1) bstep(h,dha,1)]);
        ada = min([bstep(p,dpa,1) bstep(s,dsa,1) bstep(t1,dt1a,1) bstep(t2,dt2a,1) bstep(eta,detaa,1)]);
        gap = sum(xp)+sum(vs)+sum(z1t1)+sum(z2t2)+sum(he);  mu = gap/npair;
        gap_a = sum((x+apa*dxa).*(p+ada*dpa)) + sum((v+apa*dva).*(s+ada*dsa)) ...
              + sum((z1+apa*dz1a).*(t1+ada*dt1a)) + sum((z2+apa*dz2a).*(t2+ada*dt2a)) ...
              + sum((h+apa*dha).*(eta+ada*detaa));
        sigma = (gap_a/npair/max(mu,eps))^3;

        % corrector (centering sigma*mu + second-order correction)
        r10 = sigma*mu - xp   - dxa.*dpa;
        r11 = sigma*mu - vs   - dva.*dsa;
        r12 = sigma*mu - z1t1 - dz1a.*dt1a;
        r13 = sigma*mu - z2t2 - dz2a.*dt2a;
        rHE = sigma*mu - he   - dha.*detaa;
        rc = reduce_rhs(A,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,s,v,t1,t2,z1,z2,x) ...
             + G'*((rHE + eta.*rH)./h);
        dx = cg_step(A, A2, D, Dd, G, G2, ge, rc, tol_cg);
        dw = A*dx - r9;  dz1 = dw - r7;  dz2 = r8 - dw;
        dv = r6 - dx;    ds = (r11 - s.*dv)./v;   dy = r3 - ds;
        dp = (r10 - p.*dx)./x;
        dt1 = (r12 - t1.*dz1)./z1;  dt2 = (r13 - t2.*dz2)./z2;
        dw1 = r4 + dt1;  dw2 = r5 - dt2;  dq = r2 + dw1 + dw2;
        dh = G*dx - rH;  deta = (rHE - eta.*dh)./h;

        % the smoothing removes the vertex kink, so no vertex step control is
        % needed; z1,z2 keep w inside [b1,bE].
        ap = min([bstep(x,dx,tau) bstep(v,dv,tau) bstep(z1,dz1,tau) bstep(z2,dz2,tau) bstep(h,dh,tau)]);
        ad = min([bstep(p,dp,tau) bstep(s,ds,tau) bstep(t1,dt1,tau) bstep(t2,dt2,tau) bstep(eta,deta,tau)]);

        x=x+ap*dx; w=w+ap*dw; v=v+ap*dv; z1=z1+ap*dz1; z2=z2+ap*dz2; h=h+ap*dh;
        p=p+ad*dp; s=s+ad*ds; t1=t1+ad*dt1; t2=t2+ad*dt2; eta=eta+ad*deta;
        y=y+ad*dy; q=q+ad*dq; w1=w1+ad*dw1; w2=w2+ad*dw2;

        it = it + 1;

        d = A*x;  [sup,dg] = membership_smooth(d, b, mu_min);  gradf = A'*dg;
        gap = sum(x.*p)+sum(v.*s)+sum(z1.*t1)+sum(z2.*t2)+sum(h.*eta);
        cp = max([ norm(U-x-v)/(1+norm(U)), ...
                   norm(w-A*x)/(1+norm(x)+norm(w)), ...
                   norm(bL+z1-w)/(1+norm(bL)), ...
                   norm(bE-z2-w)/(1+norm(bE)), ...
                   norm(gradf-p-y-A'*q-G'*eta)/(1+norm(gradf)), ...
                   norm(G*x-h)/(1+norm(h)), ...
                   gap/(2*sum(sup)+1) ]);

        means = struct_means(d, sizes);
        if max(abs(means - means_prev)) < tol_plan, stable = stable+1; else, stable = 0; end
        means_prev = means;

        if verbose
            fprintf('it %3d  phi=%.6g  gap=%.3g  cp=%.3e  means=%s\n', ...
                    it, sum(sup), gap, cp, mat2str(round(means,2)));
        end

        if cp < tol_kkt
            criterion = sprintf('KKT<%.0e (it %d)', tol_kkt, it);  break;
        elseif stable >= win
            criterion = sprintf('plan stable (<%.2f Gy, it %d)', tol_plan, it);  break;
        end
    end

    d = A*x;  sup = membership_smooth(d, b, mu_min);
    out.x=x; out.dose=d; out.phi=sum(sup); out.iter=it; out.time=etime(clock,t0);
    out.cp=cp; out.criterion=criterion; out.means=struct_means(d,sizes);
    out.sizes=sizes; out.names=names; out.shape=shape; out.protocol=protocol.name;
    out.dT=dT; out.Dp=Dp;
    if verbose
        fprintf('>> %s | %s | %d it | %.0f s | stop: %s\n', ...
                protocol.name, shape, it, out.time, criterion);
    end
end


% ===================== protocol-driven construction =====================
function A = full_matrix(data, pattern, exclude)
% Dose matrix of the structure matching PATTERN (excluding '(mean)' and EXCLUDE).
    nm = {data.matrix.Name};
    hit = ~cellfun('isempty', regexpi(nm, pattern, 'once'));
    if ~isempty(exclude), hit = hit & cellfun('isempty', regexpi(nm, exclude, 'once')); end
    hit = hit & cellfun('isempty', regexpi(nm, '\(mean\)', 'once'));
    idx = find(hit, 1);
    if isempty(idx), error('Structure not found: %s', pattern); end
    A = sparse(double(data.matrix(idx).A));
end

function b = build_b(P, shape, szt, szoar, dT)
% Fuzzy numbers per structure. Target: bilateral around dT. OAR: left-anchored
% (triangular peak at 0; trapezoidal plateau [0, goal]); support bE above the
% achievable dose so feasible plans stay inside the support.
    T = P.target;
    if strcmp(shape,'triangular')
        b = zeros(szt+sum(szoar), 3);
        b(:,1) = col([T.lo*dT], [0], szt, szoar);
        b(:,2) = col([T.peak*dT], [0], szt, szoar);
        b(:,3) = col([T.hi*dT], [P.oars.support], szt, szoar);
    else
        b = zeros(szt+sum(szoar), 4);
        b(:,1) = col([T.lo*dT],      [0], szt, szoar);
        b(:,2) = col([T.core_lo*dT], [0], szt, szoar);
        b(:,3) = col([T.core_hi*dT], [P.oars.goal], szt, szoar);
        b(:,4) = col([T.hi*dT],      [P.oars.support], szt, szoar);
    end
end

function [dT, Dp] = target_scale(problem, T, P)
% Target reference dose dT (the target's TROTS Objective) and prescription Dp
% (parsed from the target name, e.g. 'PTV 0-46 Gy' -> 46). Protocol fields
% dT/Dp override the database if provided.
    pnm = {problem.Name};
    keep = ~cellfun('isempty', regexpi(pnm, T.pattern, 'once'));
    if ~isempty(T.exclude), keep = keep & cellfun('isempty', regexpi(pnm, T.exclude, 'once')); end
    idx = find(keep, 1);
    if isempty(idx), error('Target not found in problem list: %s', T.pattern); end
    dT = problem(idx).Objective;
    tok = regexp(problem(idx).Name, '(\d+(?:\.\d+)?)\s*Gy', 'tokens');
    if ~isempty(tok), Dp = str2double(tok{end}{1}); else, Dp = dT; end
    if isfield(P,'dT') && ~isempty(P.dT), dT = P.dT; end
    if isfield(P,'Dp') && ~isempty(P.Dp), Dp = P.Dp; end
end

function v = col(tval, oarvals, szt, szoar)
% Stack a per-structure scalar into a voxel-wise column.
    v = tval*ones(szt,1);
    for j = 1:numel(szoar), v = [v; oarvals(min(j,end))*ones(szoar(j),1)]; end
end


% ===================== interior-point building blocks =====================
function rc = reduce_rhs(A, r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13, s,v,t1,t2,z1,z2,x)
    rw1 = r4 + (r12 + (r7+r9).*t1)./z1;
    rw2 = r5 - (r13 - (r8+r9).*t2)./z2;
    rc  = r1 + r3 + (s.*r6 - r11)./v + A'*(r2 + rw1 + rw2) + r10./x;
end

function dx = cg_step(A, A2, D, Dd, G, G2, ge, rhs, tol_cg)
% Solve (diag(Dd) + A'diag(D)A + G'diag(ge)G) dx = rhs matrix-free by
% diagonally preconditioned CG (never forms the dense normal matrix). The G
% term is skipped entirely when there are no linear constraints.
    if isempty(G)
        H = @(vv) A'*(D.*(A*vv)) + Dd.*vv;
        M = A2'*D + Dd;
    else
        H = @(vv) A'*(D.*(A*vv)) + Dd.*vv + G'*(ge.*(G*vv));
        M = A2'*D + Dd + G2'*ge;
    end
    [dx,~] = pcg(H, rhs, tol_cg, 500, @(r) r./M);
end

function G = find_constraints(data, problem, N)
% Discover PTV min/max dose-constraint matrices from the problem list:
% a PTV entry that is a constraint with zero objective. Minimise==0 is a
% minimum-dose matrix (A x >= 0); Minimise==1 is a maximum-dose matrix
% (A x <= 0, stacked as -A x >= 0). All returned as G with G x >= 0.
    G = [];
    for i = 1:numel(problem)
        pr = problem(i);
        if isempty(pr.Name) || isempty(pr.Objective), continue; end
        if ~isempty(regexpi(pr.Name,'PTV','once')) && pr.IsConstraint==1 && pr.Objective==0
            Ai = sparse(double(data.matrix(pr.dataID).A));
            if pr.Minimise==0, G = [G; Ai]; else, G = [G; -Ai]; end
        end
    end
    if isempty(G), G = sparse(0,N); end
    G = sparse(G);
end

function a = bstep(wv, dw, fac)
    a = 1;  neg = dw < 0;
    if any(neg), a = min(1, fac*min(-wv(neg)./dw(neg))); end
end

function m = struct_means(d, sizes)
    m = zeros(1, numel(sizes));  o = 0;
    for j = 1:numel(sizes), m(j) = mean(d(o+1:o+sizes(j))); o = o+sizes(j); end
end

function x0 = init_point(data)
% Least-squares warm start from the TROTS initialisation fields. Each init
% matrix has its own reference dose (dref may be a vector), so the RHS pairs
% every matrix with its dose; the regularisation block targets zero.
    ids  = data.misc.InitialiseMatrixID;
    dref = data.misc.InitialiseReferenceDose;
    idB  = data.misc.InitialiseRegularisationMatrixID;
    T = [];  rhs = [];
    for k = 1:numel(ids)
        Ak = double(data.matrix(ids(k)).A);
        T  = [T; Ak];
        dk = dref(min(k, numel(dref)));       % this matrix's reference dose
        rhs = [rhs; dk*ones(size(Ak,1),1)];
    end
    B = double(data.matrix(idB).A);
    [x0,~] = lsqr([T;B], [rhs; zeros(size(B,1),1)], 1e-6, 100);
end


% ===================== C^1 smoothed fuzzy surprise =====================
function [sup, dg, dH] = membership_smooth(d, b, um)
% Surprise S(mu(d)) with derivatives. Triangular b has 3 columns, trapezoidal
% 4. S is (1/mu-1)^2 for mu>=um and a C^1 quadratic extension below um (no
% pole, no kink). Returns value, dS/dd and d2S/dd2.
    n = size(b,2);  bL = b(:,1);
    u = zeros(size(d));  dudd = zeros(size(d));
    if n == 3
        b2 = b(:,2);  bU = b(:,3);
        L = d>bL & d<b2;  R = d>b2 & d<bU;
        u(L)=(d(L)-bL(L))./(b2(L)-bL(L));  dudd(L)= 1./(b2(L)-bL(L));
        u(R)=(bU(R)-d(R))./(bU(R)-b2(R));  dudd(R)=-1./(bU(R)-b2(R));
        u(d==b2) = 1;
    else
        b2 = b(:,2);  b3 = b(:,3);  bU = b(:,4);
        L = d>bL & d<b2;  R = d>b3 & d<bU;
        u(L)=(d(L)-bL(L))./(b2(L)-bL(L));  dudd(L)= 1./(b2(L)-bL(L));
        u(R)=(bU(R)-d(R))./(bU(R)-b3(R));  dudd(R)=-1./(bU(R)-b3(R));
        u(d>=b2 & d<=b3) = 1;
    end
    [Sv, dS, d2S] = surprise(u, um);
    sup = Sv;  dg = dS.*dudd;  dH = d2S.*(dudd.^2);
end

function [Sv, dS, d2S] = surprise(u, um)
    Sv = zeros(size(u));  dS = zeros(size(u));  d2S = zeros(size(u));
    hi = u >= um;  lo = ~hi;  uh = u(hi);
    Sv(hi)  = (uh-1).^2 ./ uh.^2;
    dS(hi)  = 2*(uh-1) ./ uh.^3;
    d2S(hi) = 2*(3-2*uh) ./ uh.^4;
    S0 = (um-1)^2/um^2;  D1 = 2*(um-1)/um^3;  kap = 2*(3-2*um)/um^4;
    ul = u(lo);
    Sv(lo)  = S0 + D1*(ul-um) + 0.5*kap*(ul-um).^2;
    dS(lo)  = D1 + kap*(ul-um);
    d2S(lo) = kap;
end