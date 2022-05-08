function DF1(x, tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    G = abs(sin(0.5*pi*t))
    H = 0.75*sin(0.5*pi*t)+1.25
    g = 1+sum((x[2:end]-G).^2)
    f1 = x[1]
    f2 = g*(1-(x(1)/g)^H)
    [f1, f2]
end

function DF1_PF(tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    x = LinRange(0, 1, 1500)
    G = 0.75*sin(0.5*pi*t)
    H = G+1.25
    f1 = x
    f2 = g*(1-(x/g).^H)
    get_PF([f1, f2], false)
end


function DF2(x, tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    G = abs(sin(0.5*pi*t))
    r = 1+floor((n-1)*G)
    tmp = setdiff(1:n,r)
    g = 1+sum((x(tmp)-G).^2)
    f1 = x[r]
    f2 = g*(1-(x(r)/g)^0.5)x
    [f1, f2]
end

function DF2_PF(tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    x = LinRange(0, 1, 1500)
    g = abs(sin(0.5*pi*t))
    f1 = x
    f2 = g*(1-(x/g).^0.5)
    get_PF([f1, f2], false)
end

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

function DF4(x, tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    a = sin(0.5*pi*t)
    b = 1+abs(cos(0.5*pi*t))
    c = max(abs(a), a+b)
    H = 1.5+a
    g = 0 # 1+sum((x[2:end]-a*(x[1]/c).^(2./[2:n])).^2)
    f1 = g*abs(x(1)-a).^H
    f2 = g*abs(x(1)-a-b).^H
    [f1, f2]
end

function DF4_PF(tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    x = LinRange(0, 1, 1500)
    a = sin(0.5*pi*t)
    b = 1+abs(cos(0.5*pi*t))
    x = linspace(a,a+b)
    H = 1.5+a
    f1 = g*abs(x-a).^H
    f2 = g*abs(x-a-b).^H
    get_PF([f1, f2], false)
end

function DF5(x, tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    G = sin(0.5*pi*t)
    w = floor(10*G)
    g = 1+sum((x[2:end]-G).^2)
    f1 = g*(x(1)+0.02*sin(w*pi*x(1)))
    f2 = g*(1-x(1)+0.02*sin(w*pi*x(1)))
    [f1, f2]
end

function DF5_PF(tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    G = sin(0.5*pi*t)
    w = floor(10*G)
    f1 = g*(x+0.02*sin(w*pi*x))
    f2 = g*(1-x+0.02*sin(w*pi*x))
    get_PF([f1, f2], false)
end

function DF6(x, tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    G = sin(0.5*pi*t)
    a = 0.2+2.8*abs(G)
    y = x[2:end]-G
    g = 1+sum((abs(G)*y.^2-10*cos(2*pi*y)+10))
    f1 = g*(x1+0.1*sin(3*pi*x1)).^a
    f2 = g*(1-x1+0.1*sin(3*pi*x1)).^a
    [f1, f2]
end

function DF6_PF(tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    G=sin(0.5*pi*t)
    a=0.2+2.8*abs(G)
    f1=g*(x+0.1*sin(3*pi*x)).^a
    f2=g*(1-x+0.1*sin(3*pi*x)).^a
    get_PF([f1, f2], false)
end

function DF7(x, tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    a=5*cos(0.5*pi*t)
    tmp=1/(1+exp(a*(x(1)-2.5)))
    g=1+sum((x[2:end]-tmp).^2)
    f1=g*(1+t)/x(1)
    f2=g*x(1)/(1+t)
    [f1, f2]
end

function DF7_PF(tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    x = LinRange(0, 1, 1500)
    f1=g*(1+t)./x
    f2=g*x/(1+t)
    get_PF([f1, f2], false)
end

function DF8(x, tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    G=sin(0.5*pi*t)
    a=2.25+2*cos(2*pi*t)
    tmp=G*sin(4*pi*x(1))/(1+abs(G))
    g=1+sum((x[2:end]-tmp).^2)
    f1=g*(x(1)+0.1*sin(3*pi*x(1)))
    f2=g*(1-x(1)+0.1*sin(3*pi*x(1))).^a
    [f1, f2]
end

function DF8_PF(tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    x = LinRange(0, 1, 1500)
    a=2.25+2*cos(2*pi*t)
    f1=g*(x+0.1*sin(3*pi*x))
    f2=g*(1-x+0.1*sin(3*pi*x)).^a
    get_PF([f1, f2], false)
end

function DF9(x, tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    N=1+floor(10*abs(sin(0.5*pi*t)))
    g=1
    for i=2:n
        tmp=x(i)-cos(4*t+x(1)+x(i-1))
        g=g+tmp.^2
    end
    f1=g*(x(1)+max(0, (0.1+0.5/N)*sin(2*N*pi*x(1))))
    f2=g*(1-x(1)+max(0, (0.1+0.5/N)*sin(2*N*pi*x(1))))
    [f1, f2]
end

function DF9_PF(tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    x = LinRange(0, 1, 1500)
    N=1+floor(10*abs(sin(0.5*pi*t)));
    f1=g*(x+max(0, (0.1+0.5/N)*sin(2*N*pi*x)));
    f2=g*(1-x+max(0, (0.1+0.5/N)*sin(2*N*pi*x)));
    get_PF([f1, f2], false)
end

function DF10(x, tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    G=sin(0.5*pi*t)
    H=2.25+2*cos(0.5*pi*t)
    tmp=sin(4*pi*(x(1)+x(2)))/(1+abs(G))
    g=1+sum((x[3:end]-tmp).^2)
    f1=g*sin(0.5*pi*x(1)).^H
    f2=g*sin(0.5*pi*x(2)).^H.*cos(0.5*pi*x(1)).^H
    f3=g*cos(0.5*pi*x(2)).^H.*cos(0.5*pi*x(1)).^H
    [f1, f2, f3]
end

function DF10_PF(tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    x1,x2=meshgrid(linspace(0,1,H))
    H=2.25+2*cos(0.5*pi*t)
    f1=g*sin(0.5*pi*x1).^H
    f2=g*sin(0.5*pi*x2).^H.*cos(0.5*pi*x1).^H
    f3=g*cos(0.5*pi*x2).^H.*cos(0.5*pi*x1).^H
    get_PF([f1, f2, f3], false)
end

function DF11(x, tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    G=abs(sin(0.5*pi*t))
    g=1+G+sum((x[3:end]-0.5*G*x(1)).^2)
    y=pi*G/6+(pi/2-pi*G/3)*x(1:2)
    f1=g*sin(y(1))
    f2=g*sin(y(2))*cos(y(1))
    f3=g*cos(y(2))*cos(y(1))
    [f1, f2]
end

function DF11_PF(tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    G=abs(sin(0.5*pi*t))
    y1=pi*G/6+(pi/2-pi*G/3)*x1
    y2=pi*G/6+(pi/2-pi*G/3)*x2
    f1=g.*sin(y1)
    f2=g.*sin(y2)*cos(y1)
    f3=g.*cos(y2)*cos(y1)
    get_PF([f1, f2, f3], false)
end

function DF12(x, tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    k=10*sin(pi*t)
    tmp1=x[3:end]-sin(t*x(1))
    tmp2=abs(sin(floor(k*(2*x(1:2)-1))*pi/2))
    g=1+sum(tmp1.^2)+prod(tmp2)
    f1=g*cos(0.5*pi*x(2)).*cos(0.5*pi*x(1))
    f2=g*sin(0.5*pi*x(2)).*cos(0.5*pi*x(1))
    f3=g*sin(0.5*pi*x(1))
    [f1, f2, f3]
end

function DF12_PF(tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    x1,x2=meshgrid(linspace(0,1,H))
    k=10*sin(pi*t)
    tmp2=abs(sin(floor(k*(2*x1-1))*pi/2).*sin(floor(k*(2*x2-1))*pi/2))
    g=1+tmp2
    f1=g.*cos(0.5*pi*x2).*cos(0.5*pi*x1)
    f2=g.*sin(0.5*pi*x2).*cos(0.5*pi*x1)
    f3=g.*sin(0.5*pi*x1)
    get_PF([f1, f2, f3], false)
end

function DF13(x, tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    G=sin(0.5*pi*t)
    p=floor(6*G)
    g=1+sum((x[3:end]-G).^2)
    f1=g*cos(0.5*pi*x(1)).^2
    f2=g*cos(0.5*pi*x(2)).^2
    f3=g*sin(0.5*pi*x(1)).^2+sin(0.5*pi*x(1)).*cos(p*pi*x(1)).^2 +
        sin(0.5*pi*x(2)).^2+sin(0.5*pi*x(2)).*cos(p*pi*x(2)).^2
    [f1, f2, f3]
end

function DF13_PF(tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    G=sin(0.5*pi*t)
    p=floor(6*G)
    f1=g.*cos(0.5*pi*x1).^2
    f2=g.*cos(0.5*pi*x2).^2
    f3=g.*sin(0.5*pi*x1).^2+sin(0.5*pi*x1).*cos(p*pi*x1).^2 +
        sin(0.5*pi*x2).^2+sin(0.5*pi*x2).*cos(p*pi*x2).^2
    get_PF([f1, f2, f3], false)
end

function DF14(x, tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    G=sin(0.5*pi*t)
    g=1+sum((x[3:end]-G).^2)
    y=0.5+G*(x(1)-0.5)
    f1=g*(1-y+0.05*sin(6*pi*y))
    f2=g*(1-x(2)+0.05*sin(6*pi*x(2)))*(y+0.05*sin(6*pi*y))
    f3=g*(x(2)+0.05*sin(6*pi*x(2)))*(y+0.05*sin(6*pi*y))
    [f1, f2, f3]
end

function DF14_PF(tau, taut, nt)
    T0 = 50
    tau_tmp = max(tau+taut-(T0+1), 0)
    t = 1/nt*floor(tau_tmp/taut)

    n = length(x)

    x1,x2=meshgrid(linspace(0,1,H))
    G=sin(0.5*pi*t)
    y=0.5+G*(x1-0.5)
    f1=g.*(1-y+0.05*sin(6*pi*y))
    f2=g.*(1-x2+0.05*sin(6*pi*x2)).*(y+0.05*sin(6*pi*y))
    f3=g.*(x2+0.05*sin(6*pi*x2)).*(y+0.05*sin(6*pi*y))
    get_PF([f1, f2, f3], false)
end
