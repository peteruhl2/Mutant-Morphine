using NLsolve

npts = 1000
morphine = range(0,200,length=npts)
fitness = range(0,1,length=npts)
escape = range(0,50,length=npts)
results = zeros(length(escape),length(fitness))

function myfun!(val,M,F,B)

    ### q-r parameters
    Mh = 100
    rc = 0.16
    rm = 0.52
    qc = 1.23e-6
    qm = 0.25
    n = 8

    eta_r(M) = (M.^n)./(Mh^n+M.^n)
    eta_q(M) = 1-eta_r(M)

    r = rc + (rm - rc).*eta_r(M)
    q = qc + (qm - qc).*eta_q(M)

    # other parameters
    lambda = 3690 #%3690;
    # F = 0.2;%0.2;
    bl = 1e-9
    bh = 1e-7
    p = 2500 #%2500
    b = 0.25 #%0.005 Vitaly: 0.01 to 0.4
    # B = 30
    dt = 0.01
    dv = 23
    di = 0.7
    dc = 0.2

    ep = 3e-5
    mu = 1
    eta = 1
    EP = ep./(mu+eta.*M)

    alp = 6.7e-5
    gamma = 1
    xi = 1
    ALP = alp./(gamma + xi.*M)

    omega_base = 15
    psi = 0.1
    omega = omega_base*exp(-psi.*M)

    val =  lambda * p * (M * eta - ep + mu) * dc *
            (((qm + dt) * bl + bh * rc) * Mh ^ (2 * n) +
            ((qc + dt) * bl + bh * rm)*
            M ^ (2 * n) + ((2 * dt + qc + qm) * bl + bh *
            (rc + rm)) * Mh ^ n * M ^ n) / dt /
            (Mh ^ n + M ^ n) / (M * eta + mu)/
            ((qm + rc + dt) * Mh ^ n + (qc + rm + dt)*
            M ^ n) / dv /(omega_base * exp(-(psi * M)) * b + di *
            dc) - (-p * lambda * (((qm + dt) * bl + bh * rc) *
            Mh ^ (2 * n) + ((qc + dt) * bl + bh * rm) *
            M ^ (2 * n) + ((2 * dt + qc + qm) * bl + bh *
            (rc + rm)) * Mh ^ n * M ^ n) * (1 + B) * (-1 + F) *
            dc / dv / dt / (Mh ^ n + M ^ n) / ((qm + rc + dt) *
            Mh ^ n + (qc + rm + dt) * M ^ n) /
            (omega_base * exp(-(psi * M)) * b + di * dc * (1 + B)))
end
