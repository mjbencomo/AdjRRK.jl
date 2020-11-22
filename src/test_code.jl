# DERIVATIVE TEST
function derv_test!(solver!,arrks::AdjRRK_struct,arrks_h::AdjRRK_struct,ts::Time_struct,rk,h0,Nref)

    solver!(arrks,ts,rk)
    @unpack u,u0 = arrks

    solver!(arrks,ts,rk;lin=true)
    @unpack u_lin,u0_lin = arrks

    errs = zeros(Nref+1)
    rate = zeros(Nref)

    h    = zeros(Nref+1)
    h[1] = h0
    for n=1:Nref
        h[n+1] = h[n]/2
    end

    for n=1:Nref+1
        arrks_h.u0 = u0 + h[n].*u0_lin
        solver!(arrks_h,ts,rk)
        uh = arrks_h.u
        errs[n] = norm((uh[:,end]-u[:,end])./h[n] - u_lin[:,end])
    end

    for n=1:Nref
        rate[n] = log2(errs[n]) - log2(errs[n+1])
    end

    return errs,rate,h
end

# INNER PRODUCT TEST
function ip_test!(solver!,arrks::AdjRRK_struct,ts::Time_struct,rk)
    solver!(arrks,ts,rk)

    arrks.u0_lin = randn(length(arrks.u0))
    solver!(arrks,ts,rk;lin=true)

    arrks.uT_adj = randn(length(arrks.u0))
    solver!(arrks,ts,rk;adj=true)

    @unpack u_lin,u_adj = arrks
    ipt = u_lin[:,1]⋅u_adj[:,1] - u_lin[:,end]⋅u_adj[:,end]
    ipt = abs(ipt)
    return ipt /= norm(u_lin)*norm(u_adj[:,end])
end

# CONVERGENCE TEST
function conv_test!(solver!,arrks::AdjRRK_struct,ts::Time_struct,rk,dt0,Nref,u_true)
    dt = zeros(Nref+1)
    dt[1] = dt0
    for n=1:Nref
        dt[n+1] = dt[n]/2
    end

    errs = zeros(Nref+1)
    for n=1:Nref+1
        ts.dt = dt[n]
        solver!(arrks,ts,rk)
        @unpack u = arrks
        errs[n] = norm(u[:,end] - u_true)
    end

    rates = zeros(Nref)
    for n=1:Nref
        rates[n] = log2(errs[n]) - log2(errs[n+1])
    end
    return errs,rates,dt
end
