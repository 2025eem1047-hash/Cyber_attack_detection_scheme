function [V_est, ang_est, h_final, H_final] = run_wls_estimation(Ybus, branch, z, R, busdata)

tol = 1e-8;         % convergence tolerance on state update norm
max_iter = 50;      % max Gauss-Newton iterations
eps_V = 1e-6;       % finite diff for V
eps_ang = 1e-6;     % finite diff for angle (rad)

nb = size(Ybus,1);
nl = length(branch);
m = length(z);
nstate = nb + nb;   % [V; ang] both of length nb

% initial guess from busdata
if size(busdata,2) >= 8
    V = busdata(:,7);
    ang = deg2rad(busdata(:,8));
else
    % fallback
    V = ones(nb,1);
    ang = zeros(nb,1);
end

% Precompute Rinv
Rinv = inv(R);

% measurement function as nested function
    function h = meas_from_state(Vloc, angloc)
        % returns [Pinj; Pflow] with same ordering as your script
        % Pinj: real(Sbus) where Sbus = V*exp(j*ang) .* conj(Ybus * V*exp(j*ang))
        Vcomplex = Vloc .* exp(1j*angloc);
        Sbus = Vcomplex .* conj(Ybus * Vcomplex);
        Pinj = real(Sbus);

        Pflow = zeros(nl,1);
        for kk = 1:nl
            f = branch(kk).f;
            t = branch(kk).t;
            y = branch(kk).y;
            bsh = branch(kk).b_sh;
            if isfield(branch(kk),'tap') && ~isempty(branch(kk).tap)
                tap = branch(kk).tap;
                if tap == 0, tap = 1; end
            else
                tap = 1;
            end

            Vf = Vloc(f) * exp(1j * angloc(f));
            Vt = Vloc(t) * exp(1j * angloc(t));

            % Current from f to t (considering tap on transformer)
            I_ft = (Vf - Vt / tap) * y + 1j * (bsh/2) * Vf;
            Pflow(kk) = real(Vf * conj(I_ft));
        end

        h = [Pinj; Pflow];
    end

% Start Gauss-Newton iterations
for iter = 1:max_iter
    % compute predicted measurements and residual
    h = meas_from_state(V, ang);
    r = z - h;

    % build Jacobian H (m x nstate) via finite differences
    H = zeros(m, nstate);

    % perturb voltage magnitudes
    for i = 1:nb
        dV = eps_V * max(1, abs(V(i)));
        Vp = V; Vm = V;
        Vp(i) = Vp(i) + dV;
        Vm(i) = Vm(i) - dV;
        hp = meas_from_state(Vp, ang);
        hm = meas_from_state(Vm, ang);
        H(:, i) = (hp - hm) / (2 * dV);
    end

    % perturb angles
    for i = 1:nb
        dth = eps_ang * max(1, abs(ang(i)));
        thp = ang; thm = ang;
        thp(i) = thp(i) + dth;
        thm(i) = thm(i) - dth;
        hp = meas_from_state(V, thp);
        hm = meas_from_state(V, thm);
        H(:, nb + i) = (hp - hm) / (2 * dth);
    end

    % Gauss-Newton update: delta = (H' R^-1 H) \ (H' R^-1 r)
    A = H' * Rinv * H;
    bvec = H' * Rinv * r;

    % regularize A if ill-conditioned
    reg = 1e-12;
    delta = (A + reg * eye(nstate)) \ bvec;

    % update states (note: delta adds because we solved H'*R^-1*(z-h) )
    V = V + delta(1:nb);
    ang = ang + delta(nb+1:end);

    % check convergence
    if norm(delta) < tol
        fprintf('WLS converged in %d iterations (||delta||=%.3e)\n', iter, norm(delta));
        break;
    end

    if iter == max_iter
        warning('WLS did not converge in %d iterations (||delta||=%.3e)', max_iter, norm(delta));
    end
end

% final outputs
V_est = V;
ang_est = ang;
h_final = meas_from_state(V_est, ang_est);
H_final = H;

end
