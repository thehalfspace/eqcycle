#################################################################
#
#   SCRIPT FOR COMPARISON OF RESULTS (MFD, CUMULATIVE SLIP, ETC)
#
#################################################################

using Plots
#...........
# Plot MFD
#...........
function MwPlot(Mw500, Mw1500, Mw3000, Mw5000, Mw7000)

    hist500 = fit(Histogram, Mw500, nbins = 10)
    hist1500 = fit(Histogram, Mw1500, nbins = 10)
    hist3000 = fit(Histogram, Mw3000, nbins = 10)
    hist5000 = fit(Histogram, Mw5000, nbins = 10)
    hist7000 = fit(Histogram, Mw7000, nbins = 10)

    # Cumulative
    cum500 = cumsum(hist500.weights[end:-1:1])[end:-1:1]
    cum1500 = cumsum(hist1500.weights[end:-1:1])[end:-1:1]
    cum3000 = cumsum(hist3000.weights[end:-1:1])[end:-1:1]
    cum5000 = cumsum(hist5000.weights[end:-1:1])[end:-1:1]
    cum7000 = cumsum(hist7000.weights[end:-1:1])[end:-1:1]

    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig.add_subplot(111)

    ax.plot(hist500.edges[1][1:end-1], cum500, ".", label="0.5 km width")
    ax.plot(hist1500.edges[1][1:end-1], cum1500, "*", label="1.5 km width")
    ax.plot(hist3000.edges[1][1:end-1], cum3000, ":", label="3.0 km width")
    ax.plot(hist5000.edges[1][1:end-1], cum5000, "^", label="5.0 km width")
    ax.plot(hist7000.edges[1][1:end-1], cum7000, "+", label="7.0 km width")
    ax.set_xlabel("Moment Magnitude (Mw)")
    ax.set_ylabel("Number of Earthquakes")
    ax.set_yscale("log")
    ax.set_title("Magnitude-frequency distribution")
    ax.legend(loc="upper right")
    show()

    figname = string(path, "mfd_comp.png")
    fig.savefig(figname, dpi = 300)
end


function sliprate_compare(Vf1, Vf2, Vf3, Vf4, t1, t2, t3, t4)

    fig = PyPlot.figure()
    ax = fig.add_subplot(111)

    ax.plot(t1, Vf1,lw=2, "r-", label="Avg. Node Spacing = 10 m")
    ax.plot(t2, Vf2,lw=2,"b-", label="Avg. Node Spacing = 15 m")
    ax.plot(t3, Vf3,lw=2, "g-", label="Avg. Node Spacing = 20 m")
    ax.plot(t4, Vf4,lw=2, "k-", label="Avg. Node Spacing = 40 m")

    ax.set_xlabel(" Time (yr)")
    ax.set_ylabel("Maximum Slip Rate on Fault (m/s)")

    ax.set_yscale("log")
    ax.legend(loc="upper left")
    show()

end

function sliprate_compare2(Vf1, Vf2, Vf3, Vf4, t1, t2, t3, t4)

    gr()

    plot(t1, Vf1, label=:"Avg. Node Spacing = 10 m")
    plot!(t2, Vf2, label=:"Avg. Node Spacing = 15 m")
    plot!(t3, Vf3, label=:"Avg. Node Spacing = 20 m")
    plot!(t4, Vf4, label=:"Avg. Node Spacing = 40 m")

    xaxis!(" Time (yr)")
    yaxis!("Maximum Slip Rate on Fault (m/s)", :log10)

end

function slip_compare(df1, df2, df3, df4, flt1, flt2, flt3, flt4)

    fig = PyPlot.figure()
    ax = fig.add_subplot(111)

    ax.plot(df1, -flt1/1e3, lw=2, "r-", label="Avg. Node Spacing = 10 m")
    ax.plot(df2, -flt2/1e3, lw=2, "b-", label="Avg. Node Spacing = 15 m")
    ax.plot(df3, -flt3/1e3, lw=2, "g-", label="Avg. Node Spacing = 20 m")
    ax.plot(df4, -flt4/1e3, lw=2, "k-", label="Avg. Node Spacing = 40 m")
    ax.set_xlabel("Coseismic Slip (m)")
    ax.set_ylabel("Depth (km)")
    #  ax.set_title("Interseismic Slip History")
    ax.set_ylim([0,24])
    ax.invert_yaxis()
    show()

end

function stress_slip(s10, s20, s40, del10, del20, del40)

    fig = PyPlot.figure()
    ax = fig.add_subplot(111)

    ax.plot(del10, s10,lw=2,"b-", label="Avg. Node Spacing = 10 m")
    ax.plot(del20, s20,lw=2, "r-", label="Avg. Node Spacing = 20 m")
    ax.plot(del40, s40,lw=2, "g-", label="Avg. Node Spacing = 40 m")

    ax.set_xlabel(" Slip (m)")
    ax.set_ylabel("On-fault Stress for one event (MPa)")

    #  ax.set_yscale("log")
    ax.legend(loc="upper right")
    show()

end
