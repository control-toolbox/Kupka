tf = 1*tf

# exponential mapping
exp(p0; saveat=[]) = φ((t0, tf), x0, p0, saveat=saveat).ode_sol

# limits
p0min = 1
p0max = 3

# initial plots
function initial_plots(p0_sol)

    times = range(t0, tf, length=5) # times for wavefronts

    # plot of the flow
    plt_flow = plot()

    p0s = range(p0min, p0max, length=10)       # range for the extremals
    for i ∈ 1:length(p0s)
        sol = exp(p0s[i])
        x = [sol(t)[1] for t ∈ sol.t]
        p = [sol(t)[2] for t ∈ sol.t]
        label = i==1 ? "extremals" : false
        plot!(plt_flow, x, p, color=:blue, label=label, z_order=:back)
    end

    # plot of wavefronts
    p0s = range(p0min, p0max, length=100)   # range to get points on the wavefronts
    xs  = zeros(length(p0s), length(times))
    ps  = zeros(length(p0s), length(times))
    for i ∈ 1:length(p0s)   # get points on the wavefronts
        sol = exp(p0s[i], saveat=times)
        xs[i, :] = [z[1] for z ∈ sol.(times)]
        ps[i, :] = [z[2] for z ∈ sol.(times)]
    end

    for j ∈ 1:length(times) # plot the wavefronts
        label = j==1 ? "flow at times" : false
        plot!(plt_flow, xs[:, j], ps[:, j], color=:green, linewidth=2, label=label, z_order=:front)
    end
#    plot!(plt_flow, xlims = (-4, 1), ylims =  (p0min, p0max))
#    plot!(plt_flow, [0, 0], [p0min, p0max], color=:black,
#            xlabel="x", ylabel="p", label="x=xf", z_order=:back)
    plot!(plt_flow, formatter = x -> format(x, precision=2))

    # plot the solution
    sol = exp(p0_sol)
    x = [sol(t)[1] for t ∈ sol.t]
    p = [sol(t)[2] for t ∈ sol.t]

    plot!(plt_flow, x, p, color=:red, linestyle=:solid, linewidth=1,
            label="extremal solution", z_order=:back)
    plot!(plt_flow, [x[end]], [p[end]], seriestype=:scatter,
            markersize=5, markerstrokewidth=0.5, color=:green, label=false, z_order=:front)

    # plot the shooting function with the solution
    plt_shoot = plot(xlabel="p₀", ylabel="y")
#    plot!(plt_shoot, xlims=(p0min, p0max), ylims=(-0.5, 0.5))

    p0s = range(p0min, p0max, length=500)
    plot!(plt_shoot, p0s, S, linewidth=2, label="x(p₀)", color=:green, z_order=:front)

#    plot!(plt_shoot, [p0min, p0max], [0, 0], color=:black, label="y=0", z_order=:back)
#    plot!(plt_shoot, [p0_sol, p0_sol], [-2, 0], color=:black,
#            label="p₀ solution", linestyle=:dash, z_order=:back)
    plot!(plt_shoot, [p0_sol], [0], seriestype=:scatter,
            markersize=5, markerstrokewidth=0.5, color=:green, label=false, z_order=:front)
    plot!(plt_shoot, formatter = x -> format(x, precision=4))

    return plt_flow, plt_shoot
end

function plot_extremal!(plt, p0, lw) # add an extremal to a plot
    sol = exp(p0)
    x = [sol(t)[1] for t ∈ sol.t]
    p = [sol(t)[2] for t ∈ sol.t]
    plot!(plt, x, p, color=:red, label=false, linewidth=lw)
    return plt
end

function plot_extremals!(plt, p0s, idx) # add some extremals to a plot
    idx = idx + 1
    lw_min = 0.5
    lw_max = 2
    N = length(p0s)
    lws = range(lw_min, lw_max, N)
    @assert 0 ≤ idx ≤ N
    for i = 1:idx
        plot_extremal!(plt, p0s[i], lws[i+N-idx])
    end
    return plt
end

T(p0, d) = S(p0) + Jₛ(p0) * d # tangent equation

function plot_iterates!(plt, p0s, idx)

    # styles
    x_style = (color=:red, seriestype=:scatter, markerstrokewidth=0.5, label="", z_order=:front)
    y_style = (color=:blue, seriestype=:scatter,
            markersize=3, markerstrokewidth=0, label="", z_order=:front)
    a_style = (color=:black, linestyle=:dash, label="")
T_style = (color=:blue, z_order=:back, label="")

    #
    idx = idx + 1
    ms_min = 1
    ms_max = 5
    N = length(p0s)
    mss = range(ms_min, ms_max, N)
    @assert 0 ≤ idx ≤ N
    for i = 1:idx
        plot!(plt, [p0s[i], p0s[i]], [0, S(p0s[i])]; a_style...)
        plot!(plt, [p0s[i]], [0]; markersize=mss[i+N-idx], x_style...)
    end

    a = p0min
    b = p0max
    if 2 ≤ idx ≤ N-1 # plot the tangent
        x = p0s[idx-1]
        plot!(plt, [x], [S(x)]; y_style...)
        plot!(plt, [a, b], [T(x, a-x), T(x, b-x)]; T_style...) # tangente
    end
    if 1 ≤ idx ≤ N
        x = p0s[idx]
        plot!(plt, [x], [S(x)]; y_style...)
    end

    return plt
end

plt_flow, plt_shoot = initial_plots(p0); # initial plots

