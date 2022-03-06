function DF3(x, tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    G = sin(0.5*π*t)
    H = g+1.5
    g = 1+sum((x[2:end]-G-x[1]^H).^2)
    f1 = x[1]
    f2 = g*(1-(x(1)/g)^H)
    [f1, f2]
end

function DF3_PS(tau, taut, nt)
    T0 = 50

    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    g = 1
    H = 50

    G = sin(0.5*pi*t)
    H = G+1.5

    x1 = LinRange(0, 1, 1500)
    x2 = G .+ x1.^H
    hcat(x1, x2)
end

function DF3_PF(tau, taut, nt)
    T0 = 50

    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    g = 1
    H = 50

    x = LinRange(0, 1, 1500)
    G = sin(0.5*pi*t)
    H = G+1.5
    f1 = x
    f2 = g*(1 .- (x/g).^H)
    get_PF([f1, f2], false)
end

function get_PF(f, nondominate)
    n = length(f)
    s = length(f[1])
    h = zeros(s, n)
    for i=1:n
        fi = reshape(f[i], s, 1)
        h[:,i] = fi
    end
    h
end

function DF3_Plot()
    anim = @animate for t in 50:10:440
        ps = MOEA.DF3_PS(t, 10, 20)
        pf = MOEA.DF3_PF(t, 10, 20)
        p = Plots.plot(title="DF3", lw=2, layout=(1,2), size=(850, 400))
        scatter!(p[1],
            ps[:,1], ps[:,2],
            seriestype = :scatter,
            markersize=2, 
            color="green",
            title = "Pareto optimal solutions (t = $(lpad(t,3,"0")))",
            legend = false,
            xlims = (0, 1),
            ylims = (-1, 2)
        )
        scatter!(p[2],
            pf[:,1], pf[:,2],
            seriestype = :scatter,
            markersize=2, 
            color="blue",
            title = "true Pareto front (t = $(lpad(t,3,"0")))",
            legend = false,
            xlims = (0, 1),
            ylims = (0, 1)
        )
    end
    gif(anim, "DF3.gif", fps = 20, loop = 0)
end

function DF3_PS_A()
    anim = @animate for t in 50:10:440
        ps = MOEA.DF3_PS(t, 10, 10)
        fig = Plots.plot(
            ps[:,1], ps[:,2],
            seriestype = :scatter,
            markersize=2, 
            color=:blue,
            title = "Time $(t)",
            legend = false,
        )
        xlims!((0, 1))
        ylims!((0, 2))
    end
    gif(anim, "DF3_PS.gif", fps = 20, loop = 0)
end

function DF3_PF_A()
    anim = @animate for t in 50:10:440
        pf = MOEA.DF3_PF(t, 10, 10)
        fig = Plots.plot(
            pf[:,1], pf[:,2],
            seriestype = :scatter,
            markersize=2, 
            color=:blue,
            title = "Time $(t)",
            legend = false,
        )
    end
    gif(anim, "DF3_PF.gif", fps = 20, loop = 0)
end


