module myDraw
using Reexport
@reexport using PyCall
@reexport using PyPlot

include("dataProcessKit.jl")
using .dataProcessKit
include("compatible.jl")

#{{ curve
#[[ curve curve_controlpoint ]]
# h, x_lim, y_lim = curve(StartPoint, ControlPoint1[, ControlPoint2], EndPoint; rescale::Bool=true, ...)
# Or ... = curve(StartPoint, EndPoint; cross=CrossPoint, t=0.5)
# ControlPoint = curve_controlpoint(StartPoint, CrossPoint, EndPoint[, t=0.5]; rescale::Bool=true, ...) #Calculate control point.
# Draw a Bezier curve. Point for inputs is a numeric tuple as (x, y).
# Outputed x_lim and y_lim are two-value-tuples indicating the boundary box of curve.
# When rescale=true, the axis will be rescaled if any part of curve is outside.
# For the following parameters, see https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.PathPatch.html
#
# ControlPoint = curve_controlpoint(StartPoint, CrossPoint, EndPoint[, t=0.5])
# Use to calculate control point for a three-point Bezier curve.
# See also: barbrace
# 15 Aug 2022

export curve, curve_controlpoint
function curve(X::Tuple{Real, Real}, Y::Tuple{Real, Real}, Z::Tuple{Real, Real}; rescale::Bool=true, kw...)
    h=gca().add_patch(matplotlib.patches.PathPatch(matplotlib.path.Path(
        [X, Y, Z],
        [matplotlib.path.Path.MOVETO; fill(matplotlib.path.Path.CURVE3, 2)]); fc="none", kw...))

    #Drawing is done. Below codes are just for calculating the boundary box.
    tx=(X[1]-Y[1])/(X[1]-2*Y[1]+Z[1])
    ty=(X[2]-Y[2])/(X[2]-2*Y[2]+Z[2])
    foo=(t, i)->begin
        if t<0
            X[i]
        elseif t>1
            Z[i]
        else
            (1-t)^2*X[i] + 2*(1-t)*t*Y[i] + t^2*Z[i]
        end
    end
    xlm=extrema((X[1], Z[1], foo(tx, 1)))
    ylm=extrema((X[2], Z[2], foo(ty, 2)))
    
    if rescale
        cxlm=xlim()
        if cxlm[1]>xlm[1] || cxlm[2]<xlm[2]
            xlim(xlm)
        end
        cylm=ylim()
        if cylm[1]>ylm[1] || cylm[2]<ylm[2]
            ylim(ylm)
        end
    end
    (h, xlm, ylm)
end
function curve(X::Tuple{Real, Real}, Y::Tuple{Real, Real}; cross::Union{Tuple{Real, Real}, Nothing}=nothing, t::Real=0.5, rescale::Bool=true, kw...)
    if isnothing(cross) #Just draw a line using patch.
        h=gca().add_patch(matplotlib.patches.PathPatch(matplotlib.path.Path(
            [X, Y],
            [matplotlib.path.Path.MOVETO, matplotlib.path.Path.LINETO]); fc="none", kw...))

        #Drawing is done. Below codes are just for calculating the boundary box.
        xlm=extrema((X[1], Y[1]))
        ylm=extrema((X[2], Y[2]))
        if rescale
            cxlm=xlim()
            if cxlm[1]>xlm[1] || cxlm[2]<xlm[2]
                xlim(xlm)
            end
            cylm=ylim()
            if cylm[1]>ylm[1] || cylm[2]<ylm[2]
                ylim(ylm)
            end
        end
        (h, xlm, ylm)
    else
        C=curve_controlpoint(X, cross, Y, t)
        curve(X, C, Y; rescale=rescale, kw...)
    end
end
function curve(X::Tuple{Real, Real}, Y::Tuple{Real, Real}, Z::Tuple{Real, Real}, W::Tuple{Real, Real}; rescale::Bool=true, kw...)
    h=gca().add_patch(matplotlib.patches.PathPatch(matplotlib.path.Path(
        [X, Y, Z, W],
        [matplotlib.path.Path.MOVETO; fill(matplotlib.path.Path.CURVE4, 3)]); fc="none", kw...))

    #Drawing is done. Below codes are just for calculating the boundary box.
    tx1=(-X[1] + 2*Y[1] - Z[1] - sqrt(W[1]*X[1] - W[1]*Y[1] - X[1]*Z[1] + Y[1]^2 - Y[1]*Z[1] + Z[1]^2))/(W[1] - X[1] + 3*Y[1] - 3*Z[1])
    tx2=(-X[1] + 2*Y[1] - Z[1] + sqrt(W[1]*X[1] - W[1]*Y[1] - X[1]*Z[1] + Y[1]^2 - Y[1]*Z[1] + Z[1]^2))/(W[1] - X[1] + 3*Y[1] - 3*Z[1])
    ty1=(-X[2] + 2*Y[2] - Z[2] - sqrt(W[2]*X[2] - W[2]*Y[2] - X[2]*Z[2] + Y[2]^2 - Y[2]*Z[2] + Z[2]^2))/(W[2] - X[2] + 3*Y[2] - 3*Z[2])
    ty2=(-X[2] + 2*Y[2] - Z[2] + sqrt(W[2]*X[2] - W[2]*Y[2] - X[2]*Z[2] + Y[2]^2 - Y[2]*Z[2] + Z[2]^2))/(W[2] - X[2] + 3*Y[2] - 3*Z[2])
    foo=(t, i)->begin
        if t<0
            X[i]
        elseif t>1
            W[i]
        else
            (1-t)^3*X[i] + 3*(1-t)^2*t*Y[i] + 3*(1-t)*t^2*Z[i]+t^3*W[i]
        end
    end
    xlm=extrema((X[1], Z[1], foo(tx1, 1), foo(tx2, 1)))
    ylm=extrema((X[2], Z[2], foo(ty1, 2), foo(ty2, 2)))
    
    if rescale
        cxlm=xlim()
        if cxlm[1]>xlm[1] || cxlm[2]<xlm[2]
            xlim(xlm)
        end
        cylm=ylim()
        if cylm[1]>ylm[1] || cylm[2]<ylm[2]
            ylim(ylm)
        end
    end
    (h, xlm, ylm)
end
curve_controlpoint(X::Tuple{Real, Real}, M::Tuple{Real, Real}, Y::Tuple{Real, Real}, t::Real=0.5)=(.-M .+ X.*(t - 1)^2 .+ Y.*t^2)./(2*t*(t - 1))
#}}

#{{ beside
#[[ beside ]]
# Xcoor = beside(X_mid_coor, BinNumber[, Index(i|vec)]; binwidth=0.8)
#Calculate the x-coordinates in order to draw a group of bars/plotboxs.
#When index is a vector with length N, function will return a N x BinNumber matrix.
#See also: pbar, pboxplot
#Xiong Jieyi, Jun 5, 2015 >Feb 6, 2016

export beside
function beside(X::AbstractVector{T1}, N::Integer, I::AbstractVector{T2}=1:N; binwidth::Real=0.8) where {T1<:Real,T2<:Integer}
    all(I.<=N) || error("Index should no larger than bin number.")
    wing=(length(X)>1 ? minimum(diff(X)) : 1)*binwidth*0.5
    XL=X-wing
    binlen=N==1 ? 0 : 2*wing/N
    hcat(map(x->binlen*(x-1)+binlen/2 .+ XL, I)...)'
end
beside(X::AbstractVector{T}, N::Integer, I::Integer) where {T<:Real}=vec(beside(X,N,[I]))
#}}

#{{ pbar
#[[ pbar ]]
#barhandle, bar_Xposition = pbar(Y; barh=false, sub=(N,I), width=0.8, from=0, err=..., ab_err=false, bar_param...)
#... = pbar(Y; grp=Y_len_group|fastgrp, colors=grp_rnum_group, ...) #color by group
#... = pbar([Y1; Y2; ...]; iscumsum=false, barh=false, colors=N x 3 matrix|N vector, errbar=Array, grp=..., ...) #[1]
#... = pbar((Y1, Y2, ...); barh=false, barh_revYaxis=true, colors=N x 3 matrix|N vector|(stack1, stack2, ...), errbar=TupleForEach, legend=vector|inp(...), ...) #[2]
#
#Draw bar plot and show it, using PyPlot.bar.
#[1] Draw as stacked bar. If input is absolute bar hights rather than number of each part, set iscumsum=true. If grp is assigned, pbar will draw the fingle-nail-polish bars.
#[2] Draw as grouped bar. Y1, Y2, ... can either be vector or matrix(draw as grouped-bars). In the case of matrix, only the top item will be colored, and errbar must be a tuple.
#from=... can be a vector or a number, indicating the bottom or left position of bar. Note that different from matplotlib.bar(), `from=' will not change the top/right position.
#
#barh=true will draw as horizental bar.
#barh_revYaxis: Whether reverse Y axis. Only work when barh=true.
#
#To draw errorbar:
# ;errbar= scalar | 2xN matrix | Tuple( vector_low, vector_high )
#     If a scalar number, len(N) array-like object, or an Nx1
#     array-like object, errorbars are drawn at +/-value relative
#     to the data.
#     If a sequence of shape 2xN, errorbars are drawn at -row1
#     and +row2 relative to the data.
#For grouped bar, errbar should be a tuple with the same size as Y.
#For stacked bar, errbar only for the highest bar.
# ab_err: whether input relative values (=false, same beheave as matplotlib.errorbar) or absolute values (=true).
#
#Tip:
#To label each bar: PyPlot.xticks(1:N,c"label1,label2,...",rotation="vertical")
#To label legend: PyPlot.legend(c"legend1,legend2,...")
#Note the legend parameter only support in pbar((...); legend=...,...) model. In this case, only the group but not the stack will be legend.
#
#See also: pcolor, boxplot, bar, binobar, barbrace, beside, halfbox
#Xiong jieyi, December 31, 2014>January 27, 2015>March 18, 2015>March 26, 2015>Aug 27, 2015>Oct 10, 2015>Jun 7, 2016>Dec 8, 2016>17 Jan 2018>28 May 2018>19 Jun 2018>6 Apr 2023

#Single bar
function pbar(Y::AbstractVector{<:Real}; barh::Bool=false, sub::Tuple{Integer,Integer}=(1,1), width::Real=0.8, from::Union{Real, AbstractVector{<:Real}}=0.0, errbar=nothing, ab_err::Bool=false, barh_revYaxis::Bool=true, grp::Union{Group,fastgrp,Nothing}=nothing, cmap::Union{AbstractString, Symbol, Nothing}=:C, colors::Union{AbstractVector, Nothing}=nothing, wargs...)
    sub[2]<=sub[1] || error("sub[2] should be no larger than sub[1].")
    W=width/sub[1]
    X=collect(1:length(Y)) .- width/2 .+ W*(sub[2]-0.5)
    barfun=if barh
        plt.barh
    else
        plt.bar
    end
    if !isnothing(errbar) && ab_err
        if isa(errbar, AbstractMatrix)
            @assert(size(errbar, 1)==2, "errbar should be a 2xN matrix.")
            errbar=(errbar.-(Y')).*[1, -1]
        elseif isa(errbar, Tuple)
            @assert(length(errbar)==2, "errbar should be a ( Nx1, Nx1 ) tuple.")
            errbar=(Y.-errbar[1], errbar[2].-Y)
        else
            @error("When ab_err=true, errbar can only be a 2xN matrix or a ( Nx1, Nx1 ) tuple.")
        end
    end
    O=if isnothing(grp)
        if !isnothing(errbar)
            wargs=Any[(barh ? (:xerr) : (:yerr), errbar), (:capsize, 3), wargs...]
        end
        barfun(X, Y.-from, W, from; wargs...)
    else
        if !isa(grp, fastgrp)
            grp=fastgrp(grp)
        end
        if isnothing(colors)
            colors=ncolors(grp.grpnum, cmap)
        end
        if isnothing(errbar)
            grploop(grp, X, Y; input_grpno=true) do i, cX, cY
                barfun(cX, cY.-from, W, from; color=colors[i], wargs...)
            end
        else
            grploop(grp, X, Y; input_grpno=true) do i, cX, cY
                cerrbar=if isa(errbar, Matrix)
                    errbar[:, [i]]
                else
                    errbar[i]
                end
                cwargs=Any[(barh ? (:xerr) : (:yerr), cerrbar), wargs...]
                barfun(cX, cY.-from, W, from; color=colors[i], capsize=3, cwargs...)
            end
        end
    end
    
    barh_revYaxis && barh && plt.gca().invert_yaxis()
    (O,X)
end

#Stacked bar
function pbar(Y::AbstractMatrix{<:Real}; cmap::Union{AbstractString, Symbol, Nothing}=:C,
              colors::Union{AbstractVector, Nothing}=nothing, #color=nothing,
              errbar=nothing, iscumsum::Bool=false, barh::Bool=false,
              barh_revYaxis::Bool=true, grp::Union{Void, fastgrp, Group}=nothing, kwargs...)
    if isa(grp, Group)
        grp=fastgrp(grp)
    end
    if grp!=nothing
        basecolorls=map(x->[x,x,x], range(1, stop=0, length=grp.grpnum+1)[2:end-1])
    end
    # if !isnothing(color)
    #     if isa(color,AbstractString)
    #         to_rgb=matplotlib.colors.colorConverter.to_rgb
    #         colors=[to_rgb(color)...]'
    #     elseif isa(color,AbstractVector)
    #         colors=color'
    #     elseif isa(color,Tuple)
    #         colors=[color...]'
    #     else
    #         colors=color
    #     end
    #     t=range(1, stop=0, length=size(Y,1)+2)[2:end-1]
    #     colors=vcat(colors,[t t t])
    # end
    if isnothing(colors)
        if isnothing(grp)
            colors=ncolors(size(Y, 1), cmap)
        else
            colors=ncolors(size(Y, 2), cmap)
        end
    end
    t=reverse(Y, dims=1)
    sY=iscumsum ? t : cumsum(t, dims=1)
    out=Vector{Any}(undef, size(sY,1))
    outX=Vector{Any}(undef, size(sY,1))
    for i=1:size(sY, 1)
        kw_bar=Pair{Symbol, Any}[]
        if i<size(sY, 1)
            push!(kw_bar, :from => vec(sY[end-i, :]))
        end
        if errbar!=nothing && i==1
            push!(kw_bar, (barh ? (:xerr) : (:yerr)) => errbar)
            push!(kw_bar, :capsize => 3)
        end
        if isnothing(grp)
            push!(kw_bar, :color => colors[mod(i-1, length(colors))+1])
        else
            push!(kw_bar, :colors => i==1 ? colors : basecolorls[i-1])
        end
        out[i], outX[i]=pbar(vec(sY[end-i+1, :]); barh=barh, barh_revYaxis=false, kw_bar..., kwargs...)
    end
    barh_revYaxis && barh && gca().invert_yaxis()
    (out, outX)
end

#Grouped bar
function pbar(Y::Tuple;
              cmap::Union{AbstractString, Symbol, Nothing}=:C,
              colors::Union{AbstractVector, Nothing}=ncolors(length(Y), cmap),
              errbar::Tuple=(), legend=nothing, barh::Bool=false, barh_revYaxis::Bool=true, wargs...)
    plt=importpkg(:PyPlot, preloaded=true)
    out=Vector{Any}(undef, length(Y))
    outX=Vector{Any}(undef, length(Y))
    for i=1:length(Y)
        if isempty(errbar)
            cwargs=wargs
        else
            # cwargs=Any[(barh ? (:xerr) : (:yerr),errbar[i]),wargs...]
            cwargs=Any[(:barh, barh), (:errbar, errbar[i]), wargs...]
        end
        if isa(colors, Tuple)
            out[i],outX[i]=pbar(Y[i];sub=(length(Y),i), barh=barh, barh_revYaxis=false,
                                colors=colors[i], cwargs...)
        else
            out[i],outX[i]=pbar(Y[i];sub=(length(Y),i), barh=barh, barh_revYaxis=false,
                                color=getrow(colors,mod(i-1,length(colors))+1),cwargs...)
        end
    end
    if legend!=nothing
        if !isa(legend,Tuple{Tuple, Base.Iterators.Pairs})
            isa(legend,AbstractVector) || error("legend should be vector or inp(...).")
            legend=inp(legend)
        end
        legend(out, legend[1]...; legend[2]...)
    end
    barh_revYaxis && barh && gca().invert_yaxis()
    (out,outX)
end
export pbar
#}}

#{{ parallelline
#[[ parallelline ]]
# parallelline([X1 X2], Y, ...; ...)
# parallelline(X, [Y1 Y2], ...; ...)
# Draw parallel lines using plt.plot. Either X or Y should be a two-column matrix. When X is a matrix, funciton draw horizental lines. When Y is a matrix, function draw vertial lines.
#See also:
#Xiong Jieyi, 16 Jun 2020

export parallelline
function parallelline(X::AbstractMatrix{Tx}, Y::AbstractVector{Ty}, arg...; kw...) where {Tx<:Real, Ty<:Real}
    # plt=importpkg(:PyPlot, preloaded=true)
    plt=PyPlot
    if !(size(X, 2)==2 && size(X, 1)==length(Y))
        error("X should be a two-column matrix with the same row number as Y.")
    end
    xx=[X fill(NaN, size(X, 1))]
    yy=repeat(Y, 1, 3)
    plt.plot(vec(xx'), vec(yy'), arg...; kw...)
end
function parallelline(X::AbstractVector{Tx}, Y::AbstractMatrix{Ty}, arg...; kw...) where {Tx<:Real, Ty<:Real}
    # plt=importpkg(:PyPlot, preloaded=true)
    plt=PyPlot
    if !(size(Y, 2)==2 && size(Y, 1)==length(X))
        error("Y should be a two-column matrix with the same row number as X.")
    end
    yy=[Y fill(NaN, size(Y, 1))]
    xx=repeat(X, 1, 3)
    plt.plot(vec(xx'), vec(yy'), arg...; kw...)
end
#}}

#{{ halfbox
#[[ halfbox ]]
# halfbox(;ax=gca(), left=true, bottom=true, right=false, top=false,
#          hide_ticks=true)
#Show or hide each side (spines and ticks) of axis.
#When hide_ticks is ture, the ticks will also be hide.
#See also: pbar, barbrace
#Feb 4, 2017 > 19 Apr 2021

export halfbox
function halfbox(;ax=gca(), left::Bool=true, bottom::Bool=true, right::Bool=false, top::Bool=false, hide_ticks::Bool=true)
    ax.spines["left"].set_visible(left)
    ax.spines["right"].set_visible(right)
    ax.spines["top"].set_visible(top)
    ax.spines["bottom"].set_visible(bottom)

    if hide_ticks
        if left && !right
            ax.yaxis.set_ticks_position("left")
        elseif right && !left
            ax.yaxis.set_ticks_position("right")
        elseif left && right
            ax.yaxis.set_ticks_position("both")
        else
            ax.yaxis.set_ticks_position("none")
        end

        if bottom && !top
            ax.xaxis.set_ticks_position("bottom")
        elseif right && !bottom
            ax.xaxis.set_ticks_position("top")
        elseif bottom && top
            ax.xaxis.set_ticks_position("both")
        else
            ax.xaxis.set_ticks_position("none")
        end
    end
    
    nothing
end
#}}

#{{ pboxplot
#[[ pboxplot ]]
# h, tip_x, tip_y = pboxplot(Vec|(Vec1, Vec2, ...); vert=true, colors=nothing|[...], showXXX=TF, YYYprops=ds(...), ticks_kw=..., wargs...) #XXX = caps, box, fliers, means; YYY = cap, box, whisker, flier, median, mean.
# (*) h, nothing, nothing = pboxplot([G1B1 G2B1 ...; G1B2 G2B2 ...; ...] ;
#     colors=c"C0, C1, ...", style=:outline|:fill, legend=..., legend_kw=(loc="best",), ...)
# h, tip_x, tip_y[, grp_id] = pboxplot(Vec; grp=Group|fastgrp, legend_kw=(...)|nothing)
# h, tip_x, tip_y, grp_id, subgrp_id = pboxplot(Vec; grp=Group|fastgrp, subgrp=Group|fastgrp, ...)
#Draw boxplot as groups.
#The legend can also be assigned using the first output: PyPlot.legend(h, c"B1,B2,...")
#tip_x and tip_y is matrics with the same size of input, contains the coordinates for the top of each bar (with 0.02x gap), which could be used for txtplot() or barbrace(). In subgrp mode, the outputed tip_x, tip_y, grp_id and subgrp_id are row-matched. (* In this syntax, tip_x|y is not supported yet. Need further work.)
# ticks_kw=... is the keyword parameters of xticks() or yticks(), depending on vert=true or false.
#See also: pbar, beside, barbrace, halfbox
#Xiong Jieyi, Jun 5, 2015>Oct 5, 2015>Jan 11, 2016>28 Jun 2018>16 Aug 2018>13 Feb 2019>8 Mar 2019>16 May 2022

export pboxplot
function pboxplot(X; grp::Union{Group,fastgrp,Void}=nothing, subgrp::Union{Group,fastgrp,Void}=nothing, vert::Bool=true, legend_kw=NamedTuple(), style::Symbol=:outline, colors::Union{Group, Void}=style==:fill ? f"C$1".((0:size(X,1)).%10) : nothing, ticks_kw=NamedTuple(), wargs...)
    # plt=importpkg(:PyPlot, preloaded=true)
    plt=PyPlot
    lb=nothing
    if grp!=nothing
        isa(X,AbstractVector) && eltype(X)<:Real || error("When grp is given, all data should be inputed as a vector.")
        if !isa(grp, fastgrp)
            grp=fastgrp(grp)
        end
        if subgrp==nothing
            X=grpfun(x->typeof(x)[x],grp,X)
            lb=grp.grpid
        else
            if !isa(subgrp, fastgrp)
                subgrp=fastgrp(subgrp)
            end
            if colors==nothing
                colors=f"C$1".((0:subgrp.grpnum).%10)
            end
            C=Array{typeof(X)}(undef, subgrp.grpnum, grp.grpnum)
            for i=1:subgrp.grpnum, j=1:grp.grpnum
                C[i, j]=X[intersect(want(subgrp, i), want(grp, j))]
            end
            h,xx,yy=pboxplot(C; vert=vert, colors=colors, style=style, wargs...)
            if legend_kw!=nothing
                plt.legend(h, subgrp.grpid; loc="best", legend_kw...)
            end
            if vert
                plt.xticks(1:grp.grpnum, grp.grpid; ticks_kw...)
            else
                plt.yticks(1:grp.grpnum, grp.grpid; ticks_kw...)
            end
            pl=rowpile()
            broadcast((x...)->(addrow!(pl, x...), nothing), xx, yy, trsp(grp.grpid), subgrp.grpid)
            return h, value(pl)...
        end
    end
    styleD=if style==:fill
        Dict(:medianprops=>ds(linewidth=2, color="k"),
             :whiskerprops=>ds(linestyle="--"),
             :boxprops=>ds(facecolor="0.5", alpha=0.75),
             :flierprops=>ds(marker="."),
             :meanprops=>ds(marker="d", markerfacecolor="#FFFF00", markeredgecolor="k", markeredgewidth=0.5))
    else
        Dict(:medianprops=>ds(linewidth=2, color="k"),
             :whiskerprops=>ds(linestyle="--"),
             :flierprops=>ds(marker="."),
             :meanprops=>ds(marker="d", markerfacecolor="#FFFF00", markeredgecolor="k", markeredgewidth=0.5))
    end
    kwD=Dict(wargs)
    for (k, v) in kwD
        if haskey(styleD, k)
            for (kk, vv) in v
                styleD[k][kk]=vv
            end
            delete!(kwD, k)
        end
    end
    
    h = plt.boxplot(X; vert=vert, patch_artist=style==:fill, styleD..., kwD...)
    if colors!=nothing
        if style==:fill
            for (i, el) in enumerate(h["boxes"])
                el.set_facecolor(getrow(colors, i))
            end
        else
            for ob in c"boxes,medians"
                for (i, el) in enumerate(h[ob])
                    el.set(color=getrow(colors, i))
                end
            end
            for ob in c"whiskers,caps"
                for (i, el) in enumerate(h[ob])
                    el.set(color=getrow(colors, ceil(Int, i/2)))
                end
            end
        end
    end
    if lb!=nothing
        if vert
            plt.xticks(1:length(lb), map(string,lb); ticks_kw...)
        else
            plt.yticks(1:length(lb), map(string,lb); ticks_kw...)
            plt.gca().invert_yaxis()
        end
    end
    cap_top_fun=vert ? :get_ydata : :get_xdata
    N=length(h["medians"])
    cap_ys=if isempty(h["caps"])
        fill(NaN, N)
    else
        # map(x->h["caps"][x*2][cap_top_fun]()[1], 1:N) #will trigger a warning.
        map(x->getproperty(h["caps"][x*2], cap_top_fun)()[1], 1:N)
    end
    t=if vert
        gap=plt.ylim()|>x->0.02*(x[2]-x[1])
        (1:N, cap_ys.+gap)
    else
        gap=plt.xlim()|>x->0.02*(x[2]-x[1])
        (cap_ys.+gap, 1:N)
    end
    if lb==nothing
        return h, t...
    else
        return h, t..., lb
    end
end

function pboxplot(X::Matrix{T}; colors::Group=f"C$1".((0:size(X,1)).%10), vert::Bool=true, style::Symbol=:outline, legend::Union{Nothing, AbstractVector}=nothing, legend_kw=NamedTuple(), ticks_kw=NamedTuple(), wargs...) where {T<:Union{AbstractVector,Any}}
    if T==Any && !all(x->isa(x, AbstractVector), X)
        error("Some elements of input matrix are not a vector.")
    end
    styleD=if style==:fill
        Dict(:medianprops=>ds(linewidth=2, color="k", solid_capstyle="butt"),
             :boxprops=>ds(alpha=0.75),
             :whiskerprops=>ds(linestyle="--"),
             :flierprops=>ds(marker="."),
             :meanprops=>ds(marker="d", markerfacecolor="#FFFF00", markeredgecolor="k", markeredgewidth=0.5))
    else
        Dict(:medianprops=>ds(linewidth=2, color="k"),
             :whiskerprops=>ds(linestyle="--"),
             :flierprops=>ds(marker="."),
             :meanprops=>ds(marker="d", markerfacecolor="#FFFF00", markeredgecolor="k", markeredgewidth=0.5))
    end
    kwD=Dict(wargs)
    for (k, v) in kwD
        if haskey(styleD, k)
            for (kk, vv) in v
                styleD[k][kk]=vv
            end
            delete!(kwD, k)
        end
    end
    
    cap_top_fun=vert ? :get_ydata : :get_xdata
    xx=beside(1:size(X,2),size(X,1))
    # plt=importpkg(:PyPlot, preloaded=true)
    plt=PyPlot
    lim=RefAny()
    bps, caps_y=mapr(1:size(X,1)) do i
        ccol=getrow(colors,i)
        bp=plt.boxplot(vec(X[i,:]); positions=vec(xx[i,:]), widths=0.6/size(X,1),
                       vert=vert, patch_artist=style==:fill, styleD..., kwD...)
        clim=vert ? plt.ylim() : plt.xlim()
        if i==1
            lim.x=clim
        else
            lim.x=(min(lim.x[1], clim[1]), max(lim.x[2], clim[2]))
        end
        if style==:fill
            for el in bp["boxes"]
                el.set_facecolor(ccol)
            end
        else
            for ob in c"boxes,medians,whiskers,caps,fliers"
                for el in bp[ob]
                    el.set(color=ccol, mec=ccol)
                end
            end
        end
        cap_ys=if isempty(bp["caps"])
            fill(NaN, size(X, 2))
        else
            # map(x->bp["caps"][x*2][cap_top_fun]()[1], 1:size(X, 2))
            map(x->getproperty(bp["caps"][x*2], cap_top_fun)()[1], 1:size(X, 2))
        end
        (bp["boxes"][1], cap_ys)
    end
    # V=vcat(X...)
    # lim=(minimum(V),maximum(V))
    if vert
        plt.ylim(lim.x)
        plt.xlim([0.5,size(X,2)+0.5])
        plt.xticks(1:size(X,2), 1:size(X,2); ticks_kw...)
    else
        plt.xlim(lim.x)
        plt.ylim([0.5,size(X,2)+0.5])
        plt.yticks(1:size(X,2), 1:size(X,2); ticks_kw...)
        plt.gca().invert_yaxis()
    end
    if legend!=nothing
        plt.legend(bps, legend; loc="best", legend_kw...)
    end
    return bps, nothing, nothing
    #Below code is not work so far. Need fix. 23 Nov 2021
    # t=if vert
    #     gap=plt.ylim()|>x->0.02*(x[2]-x[1])
    #     (1:N, cap_ys.+gap)
    # else
    #     gap=plt.xlim()|>x->0.02*(x[2]-x[1])
    #     (cap_ys.+gap, 1:N)
    # end
    # return bps, t...
end
#}}

#{{ binobar
#[[ binobar ]]
#(bar_handle,label_handle, tip_x, tip_y)=binobar(X::Vector, N::Vector; alpha=0.05, shownum=true, barh=false, text_kw=(...), text_fun=barh ? f"$1 / $2" : f"\$\\frac{$1}{$2}\$", pbar_params...)
#...=binobar(BitVector1, BitVector2, ...)
#...=binobar((X1,X2,...),(N1,N2,...);...)
#   or binobar(Vector[...],Vector[...],...) or binobar({...},{...},...)
#   #Draw as group. Xn and Nn are numeric vector.
#...=binobar((BitVector1, BitVector2, ...),(BitVector1, BitVector2, ...)) #The same as above. Tuple can also be replaced by Vector[...].
#...=binobar(X_stacked::Matrix, N::Vector|1-row-matrix; ...) #Draw stacked bar for X_stacked/N. (*)
#Draw ratio and binomial confidential intervals using pbar.
#When shownum=true, text_kw assign the parameter of PyPlot.text().
#tip_x and tip_y is matrics with the same size of input, contains the coordinates for the top of each bar (with 0.02x gap), which could be used for txtplot() or barbrace(). In subgrp mode, the outputed tip_x, tip_y, grp_id and subgrp_id are row-matched.
#(*) In stacked mode, error bars are only for the total proportions.
#See also: pbar, bar, barbrace, halfbox
#Xiong Jieyi, January 4, 2015>Aug 20, 2015>Aug 27, 2015>11 Dec 2017>17 Jan 2018>24 Apr 2020

export binobar
function binobar(X::AbstractVector{T1},N::AbstractVector{T2}; alpha::AbstractFloat=0.05, barh::Bool=false, shownum::Bool=true, text_inp::Tuple{Tuple, Any}=inp(), text_kw=text_inp[2], text_fun::Function=barh ? f"$1 / $2" : f"\$\\frac{$1}{$2}\$", wargs...) where {T1<:Integer,T2<:Integer}
    #Using below to get 95% confidential intervals.
    # @pkgfun(HypothesisTests, confint, BinomialTest)
    V=map(X, N) do x,n
        ## I replaced HypothesisTests.jl to R function for confidential interval calculation, since HypothesisTests.jl cannot calculate when x=0, but R function can. 21 Apr 2020
        # if x>0 && n>0 && x<n
        #     confint(BinomialTest(x,n), alpha)
        # else
        #     # warn("0 bin cannot estimate confidential interval.")
        #     n==0 ? (0.0, 0.0) : (x/n, x/n)
        # end
        if n>0
            Tuple(binomtest(x, n; conf!level=1.0-alpha)[2])
        else
            (0.0, 0.0)
        end
    end
    VL=map(x->x[1],V)
    VH=map(x->x[2],V)
    VM=X./N
    VM[isnan.(VM)].=0.0
    bar_handle=pbar(VM; errbar=[vec(VM-VL)';vec(VH-VM)'], barh=barh, wargs...)
    # plt=importpkg(:PyPlot, preloaded=true)
    plt=PyPlot
    if shownum
        if barh
            gap=plt.xlim()|>x->0.02*(x[2]-x[1])
            label_handle=map((x,y,d,n)->plt.text(x.+gap, y, text_fun(d, n); horizontalalignment="left",verticalalignment="center",text_kw...), VH,collect(1:length(X)),X,N)
            return (bar_handle, label_handle, VH.+5.5*gap, 1:length(X))
        else
            gap=plt.ylim()|>x->0.02*(x[2]-x[1])
            label_handle=map((x,y,d,n)->plt.text(x, y.+gap, text_fun(d, n); horizontalalignment="center",verticalalignment="bottom",text_kw...), collect(1:length(X)),VH,X,N)
            return (bar_handle, label_handle, 1:length(X), VH.+5.5*gap)
        end
    else
        t=if barh
            gap=plt.xlim()|>x->0.02*(x[2]-x[1])
            (VH.+gap, 1:length(N))
        else
            gap=plt.ylim()|>x->0.02*(x[2]-x[1])
            (1:length(N), VH.+gap)
        end
        return (bar_handle, nothing, t...)
    end
end
function binobar(cX::Tuple{Vararg{Vector{T}}},cN::Tuple{Vararg{Vector{T}}}; alpha::AbstractFloat=0.05, barh::Bool=false, shownum::Bool=true, text_inp::Tuple{Tuple, Any}=inp(), text_kw=text_inp[2], text_fun::Function=barh ? f"$1 / $2" : f"\$\\frac{$1}{$2}\$", wargs...) where {T<:Integer}
    (length(cX)==length(cN)) && all(length(cX[i])==length(cN[i]) for i=1:length(cX)) || error("Inconsistent X and N lengths.")
    I,X=vcatr_with(1:length(cX),cX...)
    N=vcat(cN...)
    V=map(X, N) do x,n
        if n>0
            Tuple(binomtest(x, n; conf!level=1.0-alpha)[2])
        else
            (0.0, 0.0)
        end
    end
    VL=map(x->x[1],V)
    VH=map(x->x[2],V)
    VM=X./N
    VM[isnan.(VM)].=0.0
    VE=[vec(VM-VL) vec(VH-VM)]
    (cVM,cVE)=grpfun((x,y)->(Any[x],Any[y']),fastgrp(I,1:length(cX)),VM,VE,multi_output=true)
    bar_handle,barx=pbar(tuple(cVM...); errbar=tuple(cVE...), barh=barh, wargs...)
    plt=PyPlot
    if shownum
        if barh
            gap=plt.xlim()|>x->0.02*(x[2]-x[1])
            label_handle=map((x,y,d,n)->plt.text(x.+gap, y, text_fun(d, n); horizontalalignment="left",verticalalignment="center",text_kw...), VH,[barx...;],X,N)
            return (bar_handle, label_handle, VH.+5.5*gap, 1:length(X))
        else
            gap=plt.ylim()|>x->0.02*(x[2]-x[1])
            label_handle=map((x,y,d,n)->plt.text(x, y.+gap, text_fun(d, n); horizontalalignment="center",verticalalignment="bottom",text_kw...), [barx...;],VH,X,N)
            return (bar_handle, label_handle, 1:length(X), VH.+5.5*gap)
        end
    else
        t=if barh
            gap=plt.xlim()|>x->0.02*(x[2]-x[1])
            (VH.+gap, 1:N)
        else
            gap=plt.ylim()|>x->0.02*(x[2]-x[1])
            (1:N, VH.+gap)
        end
        return (bar_handle, nothing, t...)
    end
end
function binobar(CV::Tuple{Vararg{Union{BitVector,Vector{Bool}}}}...;wargs...)
    X=tuple([[[sum(x) for x in BV]...] for BV in CV]...)
    N=tuple([[[length(x) for x in BV]...] for BV in CV]...)
    binobar(X,N;wargs...)
end
binobar(CV1::Vector{T},CV::Vector{T}...;wargs...) where {T<:Vector{Bool}}=binobar(map(x->tuple(x...),tuple(CV1,CV...))...;wargs...)
binobar(CV1::Vector{BitVector},CV::Vector{BitVector}...;wargs...)=binobar(map(x->tuple(x...),tuple(CV1,CV...))...;wargs...)
binobar(CV1::Vector{Any},CV::Vector{Any}...;wargs...)=binobar(map(x->tuple(x...),tuple(CV1,CV...))...;wargs...)
function binobar(BV::BitVector...;wargs...)
    X=[sum(x) for x in BV]
    N=[length(x) for x in BV]
    binobar(X,N;wargs...)
end
binobar(CV::Vector{Bool};wargs...)=binobar(bitpack(CV);wargs...)

#Stacked bar, added in 24 Apr 2020
function binobar(X2::AbstractMatrix{T1},N::AbstractVector{T2}; alpha::AbstractFloat=0.05, barh::Bool=false, shownum::Bool=true, text_inp::Tuple{Tuple, Any}=inp(), text_kw=text_inp[2], text_fun::Function=barh ? f"$1 / $2" : f"\$\\frac{$1}{$2}\$", wargs...) where {T1<:Integer,T2<:Integer}
    X=vec(sum(X2, dims=1))
    V=map(X, N) do x,n
        if n>0
            Tuple(binomtest(x, n; conf!level=1.0-alpha)[2])
        else
            (0.0, 0.0)
        end
    end
    VL=map(x->x[1],V)
    VH=map(x->x[2],V)
    VM=X./N
    VM[isnan.(VM)].=0.0
    VM2=X2./(N')
    VM2[isnan.(VM2)].=0.0
    bar_handle=pbar(VM2; errbar=[vec(VM-VL)';vec(VH-VM)'], barh=barh, wargs...)
    # plt=importpkg(:PyPlot, preloaded=true)
    plt=PyPlot
    if shownum
        if barh
            gap=plt.xlim()|>x->0.02*(x[2]-x[1])
            label_handle=map((x,y,d,n)->plt.text(x.+gap, y, text_fun(d, n); horizontalalignment="left",verticalalignment="center",text_kw...), VH,collect(1:length(X)),X,N)
            return (bar_handle, label_handle, VH.+5.5*gap, 1:length(X))
        else
            gap=plt.ylim()|>x->0.02*(x[2]-x[1])
            label_handle=map((x,y,d,n)->plt.text(x, y.+gap, text_fun(d, n); horizontalalignment="center",verticalalignment="bottom",text_kw...), collect(1:length(X)),VH,X,N)
            return (bar_handle, label_handle, 1:length(X), VH.+5.5*gap)
        end
    else
        t=if barh
            gap=plt.xlim()|>x->0.02*(x[2]-x[1])
            (VH.+gap, 1:N)
        else
            gap=plt.ylim()|>x->0.02*(x[2]-x[1])
            (1:N, VH.+gap)
        end
        return (bar_handle, nothing, t...)
    end
end
function binobar(X2::AbstractMatrix{T1}, N::AbstractMatrix{T2}; kw...) where {T1<:Integer, T2<:Integer}
    size(N, 1)==1 || error("The 2nd input N should be a vector or an one-row-matrix.")
    binobar(X2, vec(N); kw...)
end
#Considering support below syntax:
# binobar(M::Matrix{BitVector};wargs...)=binobar(splitdim(M);wargs...)
# binobar(M::Matrix{Any};wargs...)=binobar(splitdim(M);wargs...)
# binobar{T<:Vector{Bool}}(M::Matrix{T};wargs...)=binobar(splitdim(M);wargs...)

#}}

#{{ barbrace
#[[ barbrace ]]
#barbrace((X1,X2), Y|(Y1,Y2)[, Ytop], text=""; line_inp=inp("k-"), text_kw=(...), adjust_axis=true, gap=0.01, vert=false, balance="="|"<"|">", balanceP=3, direction=:auto|:up|:dowm|:left|:right)
#barbrace([X1,X2], [Y1,Y2]...) is also allowed, functional identical.
#barbrace([ "e|t|et" ]; ...): Use mouse. Need 2 or 3 clicks. Only work in X11 mode. e=equal Ys; t=assign Ytop.
#Draw bar brace.
#When vert=true, the Y means X-axis and X means Y-axis in above syntax.
#gap=... assigned the relative gap between barbrace and text.
#balance="<" or ">" will draw inbalance vertical lines, while the longer line is balanceP times than the shorter one.
#In mouse mode, "e" means equal tails; "t" means choose top at last click. Not support vert=true so far.
#See also: pairwise_statlabel, pbar, binobar, arrowann, topblank, curve
#Xiong Jieyi, January 5, 2015 > 7 Nov 2017 > 28 Feb 2019

export barbrace
function barbrace(X::Union{Tuple{Real,Real}, AbstractVector}, Y::Union{Tuple{Real,Real}, AbstractVector}, Ytop::Real, text::AbstractString=""; text_inp::Tuple=inp(), text_kw=text_inp[2], line_inp::Tuple=inp("k-"), adjust_axis::Bool=true, vert::Bool=false, gap::Real=0.01, direction::Symbol=:auto)
    if direction!=:auto
        if vert && direction!=:left && direction!=:right
            error("When vert=true, direction= can only be :left, :right or :auto.")
        elseif !vert && direction!=:up && direction!=:down
            error("When vert=false, direction= can only be :up, :down or :auto.")
        end            
    end

    # pp=importpkg(:PyPlot, preloaded=true)
    pp=PyPlot
    if vert
        line_handle=pp.plot([Y[1],Ytop,Ytop,Y[2]],[X[1],X[1],X[2],X[2]],line_inp[1]...; line_inp[2]...)
    else
        line_handle=pp.plot([X[1],X[1],X[2],X[2]],[Y[1],Ytop,Ytop,Y[2]],line_inp[1]...; line_inp[2]...)
    end

    if adjust_axis
        limfun=vert ? pp.xlim : pp.ylim
        limtop=limfun()|>x->Ytop+(x[2]-x[1])*(isempty(text) ? 0.02 : 0.07)
        if limfun()[2]<limtop
            limfun(limfun()[1], limtop)
        end
    end

    D=if direction==:up || direction==:right
        1
    elseif direction==:down || direction==:left
        -1
    else
        Ytop<0 ? -1 : 1
    end
    
    if !isempty(text)
        text_handle=if vert
            t=(pp.xlim()|>x->x[2]-x[1])*gap*D
            pp.text(Ytop+t,mean(X),text;horizontalalignment=D>=0 ? "left" : "right", verticalalignment="center", text_kw...)
        else
            t=(pp.ylim()|>x->x[2]-x[1])*gap*D
            pp.text(mean(X),Ytop+t,text;horizontalalignment="center",verticalalignment=D>=0 ? "bottom" : "top", text_kw...)
        end
        return (line_handle, text_handle)
    else
        return (line_handle, nothing)
    end
end
barbrace(X::Union{Tuple{Real,Real}, AbstractVector}, Y::Real, args...; wargs...)=barbrace(X, Y, NaN, args...; wargs...)
function barbrace(X::Union{Tuple{Real,Real}, AbstractVector}, Y::Real, Ytop::Real=NaN, args...; balance::String="=", balanceP::Real=2.0, vert::Bool=false, direction::Symbol=:auto, wargs...)
    if direction!=:auto
        if vert && direction!=:left && direction!=:right
            error("When vert=true, direction= can only be :left, :right or :auto.")
        elseif !vert && direction!=:up && direction!=:down
            error("When vert=false, direction= can only be :up, :down or :auto.")
        end            
    end
    D=if direction==:up || direction==:right
        1
    elseif direction==:down || direction==:left
        -1
    else
        Y[1]<0 ? -1 : 1
    end
    if isnan(Ytop)
        # pp=importpkg(:PyPlot, preloaded=true)
        pp=PyPlot
        limfun=vert ? pp.xlim : pp.ylim
        W=limfun()|>x->(x[2]-x[1])/50
        if D>=0
            Ytop=maximum(Y)+W
        else
            Ytop=minimum(Y)-W
        end
    else
        W=if D>=0
            Ytop-maximum(Y)
        else
            minimum(Y)-Ytop
        end
    end
    OY=if balance=="<"
        (Y, Y-W*balanceP)
    elseif balance==">"
        (Y-W*balanceP, Y)
    else
        balance=="=" || error("Unknown balance=$balance. It should be either \"=\", \">\", or \"<\".")
        (Y, Y)
    end
    barbrace(X, OY, Ytop, args...; vert=vert, direction=direction, wargs...)
end
function barbrace(X::Union{Tuple{Real,Real}, AbstractVector}, Y::Union{Tuple{Real,Real}, AbstractVector}, args...; vert::Bool=false, direction::Symbol=:auto, wargs...)
    if direction!=:auto
        if vert && direction!=:left && direction!=:right
            error("When vert=true, direction= can only be :left, :right or :auto.")
        elseif !vert && direction!=:up && direction!=:down
            error("When vert=false, direction= can only be :up, :down or :auto.")
        end            
    end
    D=if direction==:up || direction==:right
        1
    elseif direction==:down || direction==:left
        -1
    else
        Y[1]<0 ? -1 : 1
    end
    # pp=importpkg(:PyPlot, preloaded=true)
    pp=PyPlot
    limfun=vert ? pp.xlim : pp.ylim
    W=limfun()|>x->(x[2]-x[1])/50
    if D>=0
        Ytop=maximum(Y)+W
    else
        Ytop=minimum(Y)-W
    end
    barbrace(X,Y,Ytop, args...; vert=vert, wargs...)
end
function barbrace(typ::AbstractString="", args...; vert::Bool=false, wargs...)
    vert && error("In mouse-clicking mode, vert=true is not supported so far.")
    println("Click left and right points for bar brace.")
    # pp=importpkg(:PyPlot, preloaded=true)
    pp=PyPlot
    t=pp.ginput(2)
    X=(t[1][1],t[2][1])
    Y=(t[1][2],t[2][2])
    if 'e' in typ
        if Y[1]>=0
            Y=maximum(Y)|>x->(x,x)
        else
            Y=minimum(Y)|>x->(x,x)
        end
    end
    if 't' in typ
        println("Click top point for bar brace.")
        t=pp.ginput()
        Ytop=t[1][2]
        return barbrace(X, Y, Ytop, args...; wargs...)
    else
        return barbrace(X, Y, args...; wargs...)
    end
end
#}}

#{{ pairwise_statlabel
#[[ pairwise_statlabel ]]
# pairwise_statlabel([vec1, ..., vecN]; xx=1:N, ytop=maxdata, thickness=10%_axis, p_cutoff=0.05,
#                    show_ns=false, statfun=wilcoxtest, labelfun=f"p=$(1=.3g)" or "n.s.",
#                    testpair=Nx2_matrix, barbrace_kw...)
# Or: pairwise_statlabel(statfun, ...; ...)
# Label p-value in the plots using barbrace().
# xx and ytop: the x and y coordinates for barbrace ploting. Should be in the same length as data (ytop can also be a single number).
# statfun: custom function to calculate p-value, with two inputs as data vectors, and output p-value.
# thickness: the y-thickness of a label, used to avoid text overlaps.
# show_ns and p_cutoff: Whether to show non-significant comparsions.
# labelfun: custom function for the label text. Input is a p-value and output is a string.
# testpair: a N x 2 index matrix used to assign the tests. By default all combinations will be tested.
# See also: barbrace, ax_height, grpvec
# Xiong Jieyi, 2 May 2022

export pairwise_statlabel
function pairwise_statlabel(X::AbstractVector{<:Union{AbstractVector, Any}}; xx=1:length(X), 
                            thickness::Real=ax_height(0.1),
                            ytop::Union{AbstractVector{<:Real}, Real}=maximum.(X).+ax_height(0.03),
                            show_ns::Bool=false, p_cutoff=0.05,
                            labelfun::Function=p->p<p_cutoff ? f"p=$(1=.3g)"(p) : "n.s.",
                            statfun=(x, y)->wilcoxtest(x, y), 
                            testpair::Union{Nothing, AbstractMatrix{Int}}=nothing,
                            barbrace_kw...)
    if isa(ytop, Real)
        ytop=fill(ytop, length(X))
    end
    if !(length(X)==length(xx) && length(X)==length(ytop))
        error("data, xx and ytop have different lengths.")
    end
    mxi=if isnothing(testpair)
        buf=Int[]
        for i=1:length(X)-1
            for j=i+1:length(X)
                push!(buf, i)
                push!(buf, j)
            end
        end
        trsp(reshape(buf, 2, :))
    else
        if !(size(testpair, 2)==2)
            error("testpair should be N x 2 matrix.")
        end
        testpair
    end
    mxi=sort(mxi, dims=2)
    ri=sortri((abs.(mxi[:, 1].-mxi[:, 2]), ytop[mxi[:, 2]], ytop[mxi[:, 1]], median.(X)))
    mxi=mxi[ri, :]
    for (ai, bi) in eachrow(mxi)
        p=statfun(X[ai], X[bi])::Real
        if show_ns || p<=p_cutoff
            cytop=maximum(ytop[ai:bi])
            barbrace((xx[ai], xx[bi]), cytop, labelfun(p)::AbstractString, barbrace_kw...)
            for i=ai:bi
                ytop[i]=cytop+thickness
            end
        end
    end
end
pairwise_statlabel(statfun::Function, arg...; kw...)=pairwise_statlabel(arg...; statfun=statfun, kw...)
#}}

#{{ topblank
#[[ topblank ]]
# topblank(0.15)
#Add some space on the top of axis by adjusting the ylim.
#See also: barbrace, txtplotrl, ticklabelreplace
#Xiong Jieyi, 7 Nov 2017

export topblank
function topblank(P::Real=0.15)
    # pp=importpkg(:PyPlot, preloaded=true)
    pp=PyPlot
    ylm=pp.ylim()
    pp.ylim(ylm[1], ylm[2]+(ylm[2]-ylm[1])*P)
end
#}}

#{{ ticklabelreplace
#[[ ticklabelreplace ]]
# ticklabelreplace(:x|:y, "tick1"|ticknum => "newtick1", ...)
# Modify any automatic tick label to your text.
# Following the matplotlib style, in the input, hypen '-' will be transferred to math minus Char(0x2212) automatically.
# e.g. ticklabelreplace(:x, "-10.0"=>"<-10.0") or ticklabelreplace(:x, -10.0=>"<-10.0")
# If the Pair.first is a numeric and not in the ticks, it will be added to the ticks.
# See also: txtplot, txtplotrl, topblank
#Xiong Jieyi, 22 Apr 2020 > 13 Jan 2021 > 20 May 2022

export ticklabelreplace
function ticklabelreplace(ax::Symbol, prs::Pair{T1, T2}...) where {T1<:Union{Real, AbstractString}, T2<:AbstractString}
    prs=map(prs) do pr
        (isa(pr.first, AbstractString) ? replace(pr.first, '-'=>Char(0x2212)) : pr.first) => replace(pr.second, '-'=>Char(0x2212))
    end
    # plt=importpkg(:PyPlot, preloaded=true)
    plt=PyPlot
    tickf, limf=if ax==:x
        (plt.xticks, plt.xlim)
    else
        ax==:y || error("The first input should be either :x or :y .")
        (plt.yticks, plt.ylim)
    end
    plt.draw()
    axlim=limf()
    (V, t)=tickf()
    Vt=map(x->x.get_text(), t)
    # T=replace(map(x->x.get_text(), t), prs...)
    for pr in prs
        if isa(pr.first, AbstractString)
            replace!(Vt, pr)
        elseif any((l=V.==pr.first;))
            Vt[l].=pr.second
        else
            push!(V, pr.first)
            push!(Vt, pr.second)
        end
    end
    tickf(V, Vt)
    limf(axlim)
end
#}}

#{{ arrowann
#[[ arrowann ]]
#arrowann("text" => xypoint::(x, y); angle=0(0~360), linelength=10(0~100);
#         gap=1(0~100), x_lim=(low, high), y_lim=(low, high),
#         color="k", arrowstyle="-|>", arrowprops=(kw=val, ...), other_params...)
#
#Old syntax: (using absolute coordinates, no gap.)
#arrowann("text", xytext::(x,y), xypoint::(x,y); color=..., arrowstyle=..., arrowprops=..., other_params...)
#
#Draw a arrow with a text label using matplotlib.annotate.
#`gap' is the relative distance (0~100) between arrow head and the point.
#`linelength' and `gap' are 0~100 precentage of axis length. E.g. 100 means the whole length of axis.
#arrowstyle can be: -, ->, -[, -|>, <-, <->, <|-, <|-|>, ]-, ]-[, fancy, simple, wedge, |-|
#See also: barbrance, txtplot, txtplotrl
#Feb 7, 2017 > 31 Mar 2022

export arrowann
function arrowann(label::Union{AbstractString, Symbol}, xytext::Tuple, xy::Tuple,args...;
                  color::Union{AbstractString, Symbol}="k",
                  arrowstyle::Union{AbstractString, Symbol}="-|>",
                  arrowprops=NamedTuple(), wargs...)
    annotate(label, args...; xy=xy, xycoords="data",
                                xytext=xytext, textcoords="data", va="center", ha="center",
                                color=color, arrowprops=ds(color=color, arrowstyle=arrowstyle,
                                                         arrowprops...),
                                wargs...)
end
function arrowann((txt, (xx, yy))::Pair{<:Union{AbstractString, Symbol}, <:Tuple{<:Real, <:Real}},
                  ang::Real=0, dis::Real=10; 
                  gap::Real=1, x_lim::Union{<:Tuple{<:Real, <:Real}, Nothing}=nothing, 
                  y_lim::Union{<:Tuple{<:Real, <:Real}, Nothing}=nothing,
                  color::Union{AbstractString, Symbol}="k",
                  arrowstyle::Union{AbstractString, Symbol}="-|>",
                  arrowprops=NamedTuple(), wargs...)
    ang=mod(ang, 360)
    rad=ang*pi/180
    # plt=importpkg(:PyPlot, preloaded=true)
    plt=PyPlot
    xlm=@something(x_lim, plt.xlim())
    xr=(xx-xlm[1])/(xlm[2]-xlm[1])
    ylm=@something(y_lim, plt.ylim())
    yr=(yy-ylm[1])/(ylm[2]-ylm[1])
    
    dx=cos(rad)
    dy=sin(rad)
    
    xx1=(xr+gap*dx/100)*(xlm[2]-xlm[1])+xlm[1]
    yy1=(yr+gap*dy/100)*(ylm[2]-ylm[1])+ylm[1]
    
    xx2=(xr+(dis+gap)*dx/100)*(xlm[2]-xlm[1])+xlm[1]
    yy2=(yr+(dis+gap)*dy/100)*(ylm[2]-ylm[1])+ylm[1]
    
    ha= ang<45 || ang>315 ? "left" : 135<ang<225 ? "right" : "center"
    va= 45<=ang<=135 ? "bottom" : 225<=ang<=315 ? "top" : "center"
    plt.annotate(txt; xy=(xx1, yy1), xycoords="data",
                 xytext=(xx2, yy2), textcoords="data", va=va, ha=ha,
                 color=color, arrowprops=ds(color=color, arrowstyle=arrowstyle, arrowprops...),
                 wargs...)
end
#}}

#{{ pdensity
#[[ pdensity ]]
#h=pdensity(Y::Vector, PyPlot_plot_params...; grp=group_of_Y, ...)
#h=pdensity(Ys::iteration|[Y1 Y2 ...]::Matrix, PyPlot_plot_params...;
#    keepratio=false|times=[t1, t2, ...], rank=false, kde_kw=(bandwidth=..., boundary=...),
#    lines=["b-","r:",...], colors=N x 3 matrix|N vector, fillarea=false,
#    legend="1":"N"|nothing, legend_kw=(...,))
#Draw density curve using PyPlot. If times is given, each density curve will multiply each times before ploting. If keepratio is true, times will defined as length/mean_length.
#NaN and Inf in the data will be ignored with a warning.
#If rank is true, input will be converted as tiedrank/maximum firstly.
#If fillarea=true, function will fill the whole area below density curve.
#To repress legend, set legend=nothing. You can use legend(h, c"...,...") later.
#kde_kw: parameters of KernelDensity.kde. See below:
#  bandwidth: the bandwidth of the kernel. Default is to use Silverman's rule.
#  boundary: the lower and upper limits of the kde as a tuple. Due to the fourier transforms used internally, there should be sufficient spacing to prevent wrap-around at the boundaries.
#See also: kdescatter
#Xiong Jieyi, March 9, 2015 >Aug 20, 2015 >Sep 14, 2015>Sep 20, 2015>Feb 27, 2016>Mar 3, 2016>May 9, 2016>Jun 13, 2016>2 Mar 2019>23 Mar 2020

export pdensity
function __filternan__(x)
    l=isnan.(x) .| isinf.(x)
    any(l) && @warn(f"$1 of $2 data are NaNs/Infs, they are ignored."(sum(l),length(x)))
    x[.!l]
end
function pdensity(Y::Vector,args...;grp::Union{Group,fastgrp}=[], fillarea::Bool=false, kde_kw=NamedTuple(), wargs...)
    if isempty(grp)
        @pkgfun(KernelDensity, kde)
        @pkgfun(PyPlot, plot, legend, ylabel, fill=>pltfill)
        plotfun=fillarea ? pltfill : plot
        
        l=isnan.(Y) .| isinf.(Y)
        any(l) && @warn(f"$1 of $2 data are NaNs/Infs, they are ignored."(sum(l),length(Y)))
        Y=Y[.!l]
        grp=isempty(grp) ? grp : grp[.!l]
        
        k=kde(Y; kde_kw...)
        plotfun(k.x, k.density, args...;wargs...)
        ylabel("Density")
    else
        if !isa(grp,fastgrp)
            grp=fastgrp(grp)
        end
        C=grpvec(grp,Y)
        pdensity(tuple(C...),args...;legend=rowfun(string,grp.grpid),wargs...)
    end
end
function pdensity(Ys,args...;lines::Union{Vector,Tuple}=(),colors::Group=(),legend::Union{Void,Vector}=num2str(1:length(Ys)),times::Vector=ones(AbstractFloat,length(Ys)),keepratio::Bool=false,rank::Bool=false, legend_kw=NamedTuple(), fillarea::Bool=false, kde_kw=NamedTuple(), wargs...)
    Ys=map(__filternan__,Ys)
    if rank
        Ys=vcatdo(x->tiedrank(x)./length(x),Ys...)
    end
    if keepratio
        times=[length(x) for x in Ys]|>x->x./mean(x)
    end
    @pkgfun(KernelDensity, kde)
    @pkgfun(PyPlot, plot, legend=>pltlegend, ylabel, fill=>pltfill)
    plotfun=fillarea ? pltfill : plot
    h=map(1:length(Ys)) do i
        cline=isempty(lines) ? () : (lines[i],)
        ccol=isempty(colors) ? () : ((:color,
                                  getrow(colors,i)|>x->isa(x,Matrix) ? vec(x) : x),)
        k=kde(Ys[i]; kde_kw...)
        plotfun(k.x,k.density*times[i],cline...,args...;ccol...,wargs...)[1]
    end
    ylabel("Density")
    legend!=nothing && legend_kw!=nothing && pltlegend(h, legend; legend_kw...)
    h
end
pdensity(Ys::Matrix,args...;wargs...)=pdensity(tuple(eachcol(Ys)...),args...;wargs...)
#}}

#{{ heatmap

#[[ heatmap ]]
#heatmap(M::Matrix, args...; xlabel=..., ylabel=...,
#       ycluster=false, row_max_norm=false, pixel=+oo)
#Draw heatmap matrix using PyPlot.pcolor(default). Function will choose one automatically if no assignment.
#The NaN in the data will be ignored.
#If ycluster is true, M row of M will be rearranged by hierarchy cluster.
#ycluster could also be a AbstractVector, which shows the column used to cluster.
#If row_max_norm is true, M will be transfer to op(./,M,maximum(M,2)).
#When the row or column number is larger than pixel, matrix will be compressed.
#See also: boxplot, bar, Rheatmap
#Xiong Jieyi, 30 Sep 2014>8 Oct 2014

function needCluster()
    # isdefined(:PyCall)||require("PyCall")
    isdefined(:py_fastcluster)||@eval Main.PyCall.@pyimport fastcluster as py_fastcluster
    isdefined(:py_scipy_cluster_hierarchy)||@eval Main.PyCall.@pyimport scipy.cluster.hierarchy as py_scipy_cluster_hierarchy
end
export heatmap
function heatmap(M::Matrix{T},args1...;
                 row_max_norm::Bool=false,
                 pixel::Int=0,
                 ycluster::Union{Bool,AbstractVector}=false,
                 xlabel::AbstractVector=[],
                 ylabel::AbstractVector=[],
                 args2...) where {T<:Real}
    if !isempty(xlabel) && typein(xlabel)<:Real
        xlabel=num2str(xlabel)
    end
    if !isempty(ylabel) && typein(ylabel)<:Real
        ylabel=num2str(ylabel)
    end
    if row_max_norm
        # M=op(./,M,maximum(M,2))
        M=M./maxh(M)
    end
    if isa(ycluster,Bool)
        yclusterl=1:size(M,2)
    else
        yclusterl=deepcopy(ycluster)
        ycluster=true
    end
    if ycluster
        needCluster()
        l=py_scipy_cluster_hierarchy.leaves_list(py_fastcluster.average(M[:,yclusterl]))+1
        M=M[l,:]
        if !isempty(ylabel)
            ylabel=ylabel[l,:]
        end
    end
    if pixel>0
        if size(M,1)>pixel
            tM=zeros(pixel,size(M,2))
            # rang=linspace(1,size(M,1),pixel+1)
            rang=range(1, stop=size(M,1), length=pixel+1)
            for i=1:pixel
                tM[i,:]=mean(M[ceil(rang[i]):floor(rang[i+1]),:],1)
            end
            M=tM
            ylabel=rang[1:end-1]+(rang[2]-rang[1])/2
        end
        if size(M,2)>pixel
            tM=zeros(size(M,1),pixel)
            # rang=linspace(1,size(M,2),pixel+1)
            rang=range(1, stop=size(M,2), length=pixel+1)
            for i=1:pixel
                tM[:,i]=mean(M[:,ceil(rang[i]):floor(rang[i+1])],2)
            end
            M=tM
            xlabel=rang[1:end-1]+(rang[2]-rang[1])/2
        end
    end
    
    #Flip the y-axis only ylabel is given.
    if !isempty(ylabel)
        M=flipud(M)
        ylabel=flipud(ylabel)
    else
        ylabel=1:size(M,1)
    end

    PyPlot.pcolor(M)
    PyPlot.colorbar()
    PyPlot.xlim(0,size(M,2))
    PyPlot.ylim(0,size(M,1))
    if !isempty(xlabel) && length(xlabel)<=50
        PyPlot.xticks(0.5:(length(xlabel)-0.5),xlabel)
    end
    if !isempty(ylabel) && length(ylabel)<=50
        PyPlot.yticks(0.5:(length(ylabel)-0.5),ylabel)
    end
end
#}}

#{{ txtplot txtplotrl
#[[ txtplot txtplotrl ]]
# h_vec = txtplot(X_vec|X, Y_vec|Y, S_vec,...; colors=vector|matrix, fix=false, repel=false, repel_kw=(...), text_args...) #Vector input. For simplicity, one of (X_vec, Y_vec) could be a scalar (will be expand automatially).
# h = txtplot(X,Y,S,...; color='r'|[1 0 0];...) #One input. Don't use 'colors=...' here.
# ... = txtplotrl(relative_X, relative_Y, ...; rl=:xy|:x|:y) #Use relative coordinates (0-1).
# ... = txtplotrl("subtitle",...) #Draw subtitle in rx=0.5, ry=0.9 position.
# Draw PyPlot.text(X[i],Y[i],S[i]; color=colors[i] or colors[i,:]) for every text S. 
# relative_X and relative_Y should be 0~1 and it will be automatically transferred to absolute coordinates according to PyPlot.xlim() and PyPlot.ylim(). If only one of the values is relative value, set rl=:x or rl=:y.
# fix=true will prevent function to auto-expand the axes.
# repel=true will use python package adjustText to adjust the labels. repel_kw are the parameters of adjustText. See https://adjusttext.readthedocs.io/en/latest/
#   e.g. repel_kw=(arrowprops=ds(arrowstyle="-", color="k", lw=0.5)
#
# Some useful parameters for matplotlib.text (can be attached behind)
#  ha=[ ‘center’ | ‘right’ | ‘left’ ]
#  va=[ ‘center’ | ‘top’ | ‘bottom’ | ‘baseline’ ]
#  ma=[‘left’ | ‘right’ | ‘center’ ]
#  size, rotation, color, alpha
#  (see https://matplotlib.org/api/text_api.html#matplotlib.text.Text for more)
#See also: pbar, pdensity, arrowann, topblank
#Xiong Jieyi, March 9, 2015>Jun 5, 2015>Feb 3, 2016>Feb 16, 2016>Dec 9, 2016>13 Sep 2018>17 Dec 2019

export txtplot, txtplotrl
function txtplot(X::AbstractVector,Y::AbstractVector,S::Vector{T},args...; fix::Bool=false, colors::Array=[], repel::Bool=false, repel_kw::NamedTuple=NamedTuple(), wargs...) where {T<:AbstractString}
    X, Y, S=valid(X, Y, S, follower=3)
    # plt=importpkg(:PyPlot, preloaded=false)
    plt=PyPlot
    @il isnewplot=!plt.gca().has_data()
    h=Vector{Any}(undef, length(X))
    for i=1:length(X)
        ccol=if isempty(colors)
            ()
        elseif isa(colors, AbstractMatrix)
            (color=tuple(colors[i, :]...),)
        elseif isa(colors, AbstractVector) && (eltype(colors)<:Union{AbstractString, Tuple, AbstractArray} || eltype(colors)==Type{Any})
            (color=colors[i],)
        else
            (color=colors,)
        end
        h[i]=@il(plt.text(X[i],Y[i],S[i], args...;
             horizontalalignment="center",verticalalignment="center", ccol..., wargs...))
    end

    if !isempty(X) && (isnewplot || !fix)
        if isnewplot
            xmargin=(maximum(X)-minimum(X))*.05|>x->x>0 ? x : 0.5
            ymargin=(maximum(Y)-minimum(Y))*.05|>x->x>0 ? x : 0.5
            xbound=(minimum(X)-xmargin,maximum(X)+xmargin)
            ybound=(minimum(Y)-ymargin,maximum(Y)+ymargin)
        else
            @il oxbound=plt.xlim()
            xmargin=(maximum([X;oxbound[2]])-minimum([X;oxbound[1]]))*.05
            @il oybound=plt.ylim()
            ymargin=(maximum([Y;oybound[2]])-minimum([Y;oybound[1]]))*.05
            xbound=(min(minimum(X)-xmargin,oxbound[1]), max(maximum(X)+xmargin,oxbound[2]))
            ybound=(min(minimum(Y)-ymargin,oybound[1]), max(maximum(Y)+ymargin,oybound[2]))
        end
        @il plt.xlim(xbound)
        @il plt.ylim(ybound)
    end
    
    outh=[h...]
    if repel
        importpy("adjustText").adjust_text(outh; repel_kw...)
    end
    outh
end
txtplot(X::AbstractVector,Y::AbstractVector,S::AbstractVector,args...;wargs...)=txtplot(X,Y,map(string,S),args...;wargs...)
txtplot(X::Real,Y::AbstractVector,args...;wargs...)=txtplot(fill(X,length(Y)),Y,args...;wargs...)
txtplot(X::AbstractVector,Y::Real,args...;wargs...)=txtplot(X,fill(Y,length(X)),args...;wargs...)
function txtplot(X::Real,Y::Real,S::AbstractString,args...; fix::Bool=false, wargs...)
    # plt=importpkg(:PyPlot, preloaded=false)
    plt=PyPlot
    @il isnewplot=!invokelatest(plt.gca).has_data()
    @il h=plt.text(X,Y,S,args...;horizontalalignment="center",verticalalignment="center",wargs...)

    if isnewplot || !fix
        if isnewplot
            xmargin=0.5
            ymargin=0.5
            xbound=(X-xmargin,X+xmargin)
            ybound=(Y-ymargin,Y+ymargin)
        else
            @il oxbound=plt.xlim()
            xmargin=(maximum([X;oxbound[2]])-minimum([X;oxbound[1]]))*.05
            @il oybound=plt.ylim()
            ymargin=(maximum([Y;oybound[2]])-minimum([Y;oybound[1]]))*.05
            xbound=(min(X-xmargin,oxbound[1]), max(X+xmargin,oxbound[2]))
            ybound=(min(Y-ymargin,oybound[1]), max(Y+ymargin,oybound[2]))
        end
        @il plt.xlim(xbound)
        @il plt.ylim(ybound)
    end
    
    h
end
function txtplotrl(rX,rY,args...; rl::Symbol=:xy, wargs...)
    if !in(rl, [:xy, :x, :y])
        error("Invalid rl=:$rl. rl parameter should be :xy, :x or :y .")
    end
    # pp=importpkg(:PyPlot, preloaded=false)
    pp=PyPlot
    xx=if rl==:xy || rl==:x
        @il x1,x2=pp.xlim()
        if @il(pp.plt.get(pp.gca(),:xscale))=="log"
            10^(rX*(log10(x2)-log10(x1))+log10(x1))
        else
            (x2-x1)*rX+x1
        end
    else
        rX
    end
    yy=if rl==:xy || rl==:y
        @il y1,y2=pp.ylim()
        if @il(pp.plt.get(pp.gca(),:yscale))=="log"
            10^(rY*(log10(y2)-log10(y1))+log10(y1))
        else
            (y2-y1)*rY+y1
        end
    else
        rY
    end
    
    txtplot(xx,yy,args...;wargs...)
end
txtplotrl(S::AbstractString,args...;wargs...)=txtplotrl(0.5,0.9,S,args...;wargs...)

#}}

#{{ grpplot
#[[ grpplot ]]
# hs, grpid=grpplot(grp[ => order ]|fastgrp, X, Y, lines, alphas, ...;
#                   grp_kw=(kw=[val, missing, ...],...), 
#                   showlegend=true, legend_label_fun=..., legend_kw=(...), ...)
# hs=grpplot(X,Y,fastgrp,lines,alphas,...; ...)
# hs, grpid, grp2id = grpplot(grp[ => order ], X, Y, ...; grp2=grp2id[ => order]|fastgrp, grp2_kw=(mfc|mec=ncolors(N),), legend_label_fun[1|2]=..., ...) #Support two groups.
# Draw PyPlot.plot by group. lines and alphas could be a value, a vector for each group, or nothing for default.
# grp_kw and grp2_kw can be a vector or a matrix. For vector, the `nothing' element means this keyword parameter will use the default value (nothing => `None' in python). For matrix, each row will be assigned to each group as a vector. For `missing' elements, this keyword argument will be ignored.
# legend_label_fun=x->rowfun(string, x) give a chance to replace legend labels. The input is a Group and output should be a string vector with the same length. In grp2 mode, legend_label_fun1 and legend_label_fun2 can also be used to change grp1 and grp2 legend label seperately.
# See also: eachgrp, ncolors
#Xiong Jieyi, Mar 3, 2016>Mar 11, 2016>20 Dec 2019>16 Nov 2020>1 Mar 2022

export grpplot
function grpplot(G::fastgrp, X::AbstractVector, Y::AbstractVector, lines=nothing, alphas=nothing, args...; legend=G.grpid, legend_loc="best", legend_inp=((),Any[]), legend_kw=legend_inp[2], grp_inp=inp(), grp_kw=grp_inp[2], grp2=nothing, legend_label_fun::Function=x->rowfun(string, x), showlegend::Bool=true, wargs...) #Some undocumented parameters were used to be compatiable with old codes.
    if !isnothing(grp2)
        error("In grp2 mode, the first input as fastgrp is not supported so far. Consider using group => order instead.")
    end
    # pp=importpkg(:PyPlot)
    pp=PyPlot
    hs=grpfun(G,X,Y,input_grpno=true) do i,xx,yy
        if lines==nothing
            cline=()
        elseif isa(lines,AbstractVector)
            cline=(lines[i],)
        else
            cline=(lines,)
        end
        if alphas==nothing
            calpha=()
        elseif isa(alphas,AbstractVector)
            calpha=((:alpha,alphas[i]),)
        else
            calpha=((:alpha,alphas),)
        end
        cwargs=Pair{Symbol, Any}[]
        for (kw,val) in pairs(grp_kw)
            cval=if isa(val,AbstractMatrix)
                vec(val[i,:])
            else
                val[i]
            end
            if !ismissing(cval)
                push!(cwargs, kw => cval)
            end
        end
        pp.plot(xx,yy,cline...,args...;calpha..., cwargs..., wargs...)
    end
    if showlegend && legend!=nothing && legend_kw!=nothing
        pp.legend(hs,legend_label_fun(legend);loc=legend_loc,legend_kw...)
    end
    hs
end
function grpplot(G::Union{Group, Pair}, X::AbstractVector, Y::AbstractVector, args...; grp2=nothing, wargs...)
    if !isnothing(grp2)
        return _grpplot2(G, grp2, X, Y, args...; wargs...)
    end
    O=if isa(G, Pair)
        fastgrp(G.first, G.second)
    else
        fastgrp(G)
    end
    (grpplot(O,X,Y,args...;wargs...), O.grpid)
end
function _grpplot2(G::Union{Group, Pair}, G2::Union{Group, fastgrp, Pair}, X::AbstractVector, Y::AbstractVector, arg...; grp2_kw::NamedTuple=NamedTuple(), legend_kw::NamedTuple=NamedTuple(), legend_label_fun::Function=x->rowfun(string, x), legend_label_fun1::Function=legend_label_fun, legend_label_fun2::Function=legend_label_fun, showlegend::Bool=true, grp_inp=inp(), grp_kw=grp_inp[2], kw...)
    if isa(G2, Pair)
        G2=fastgrp(G2.first, G2.second)
    elseif !isa(G2, fastgrp)
        G2=fastgrp(G2)
    end
    if isa(G, Pair)
        (G, uG)=(G.first, G.second)
    else
        uG=unival(G)
    end
    if isempty(grp2_kw)
        t=ncolors(G2.grpnum)
        grp2_kw=(mfc=t, mec=t)
    end
    h=map(eachgrp(G2, G, X, Y, no=true)) do (gi, cG, cX, cY)
        cwargs=Pair{Symbol, Any}[]
        for (kw,val) in pairs(grp2_kw)
            cval=if isa(val,AbstractMatrix)
                vec(val[gi,:])
            else
                val[gi]
            end
            if !ismissing(cval)
                push!(cwargs, kw => cval)
            end
        end
        grpplot(fastgrp(cG, uG), cX, cY, arg...; showlegend=false, grp_inp=grp_inp, grp_kw=grp_kw, cwargs..., kw...)
    end
    
    #Prepare pseudolegend
    if showlegend
        M, _=if length(arg)<1 || !(isa(arg[1], AbstractString) || isa(arg[1], AbstractVector{<:AbstractString}))
            fill("-", rownum(uG)), arg
        elseif isa(arg[1], AbstractString)
            fill(arg[1], rownum(uG)), arg[2:end]
        else
            arg[1], arg[2:end]
        end
        t=filter(in(Set([:color, :mfc, :mec])), keys(grp2_kw))
        Ccolor=if isempty(t)
            fill(missing, G2.grpnum)
        else
            grp2_kw[t[1]]
        end
        nkw=[k=>[fill(v, rownum(uG)); fill(missing, G2.grpnum)] for (k, v) in kw]
        freelegend([M; fill("[]", G2.grpnum)], [legend_label_fun1(uG); legend_label_fun2(G2.grpid)];
                   each_kw=(color=[fill("0.5", rownum(uG)); Ccolor], nkw...,), legend_kw...)
    end
    (h, uG, G2.grpid)
end
#}}

#{{ pltsave
#[[ pltsave ]]
#pltsave(filename = "~/.tempfig.pdf"; format="pdf", rasterized=false, ...)
#Save a PyPlot figure as PDF file (or other formats like png, svg, ...), with default parameters fonttype=42, bbox_inches="tight" and dpi=300.
#rasterized: rasterize every axes of the current figure before saving. The rasterized status will be changed back after saving.
#See https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.savefig.html for more parameters.
#Seaborn font problem: In the case some fonts cannot be showed in pdf file by pltsave, try importpy("seaborn").set(font="Bitstream Vera Sans") before drawing.
#See also:rfig, PDFfigs, showimfile
#Xiong Jieyi, May 13, 2015 > 17 Jul 2019 > 15 Dec 2021

export pltsave
function pltsave(figfile::AbstractString=f"$1/.tempfig.pdf"(ENV["HOME"]); format::String="pdf", rasterized::Bool=false, kw...)
    if !endswith(figfile, ".$format")
        figfile=figfile*".$format"
    end
    # plt=importpkg(:PyPlot, preloaded=true)
    plt=PyPlot
    if format=="pdf"
        plt.rc("pdf", fonttype=42)
    end
    if rasterized
        axesvec=plt.gcf().get_axes()
        axis_rasterized_stat=Vector{Union{Bool, Nothing}}(undef, length(axesvec))
        for (i, ca) in enumerate(axesvec)
            axis_rasterized_stat[i]=ca.get_rasterized()
            ca.set_rasterized(true)
        end
    end
    h=plt.savefig(figfile; bbox_inches="tight", format=format, dpi=300, kw...)
    if rasterized
        for (ca, st) in zip(axesvec, axis_rasterized_stat)
            ca.set_rasterized(st)
        end
    end
    return h
end
#}}

#{{ colormap ncolors grp2color
#[[ colormap ncolors grp2color ]]
# hex_string_vec = colormap( vector[, :rainbow]; min_=..., max_=..., rev=false, ends_gap=false)
# - Convert numeric vector to color codes. Note the keywords are min_= and max_=.
# hex_string_vec = ncolors( N[, :rainbow|:C(if N<10)]; rev=false, loop=false, ends_gap=false)
# -  Output a set of N color codes. :C means C0, C1, ... colors in hex code.
# Color4Group = grp2color(Group|fastgrp(G, uG)[, :rainbow|:C(if N<10)]; rev=false, loop=false, ends_gap=false)
# Vector_of_vector | c"C0, C1, ..." = colormapv|grp2colorv(...)
# -  Map different color labels to each group elements. The output is a N x 3 matrix.
#For available cmaps, See http://matplotlib.org/examples/color/colormaps_reference.html
#:Accent :Blues  :BrBG   :BuGn   :BuPu   :CMRmap :ColormapRegistry       :Dark2  :GnBu   :Greens :Greys:Mapping :MutableMapping :OrRd   :Oranges        :PRGn   :Paired :Pastel1        :Pastel2        :PiYG :PuBu    :PuBuGn :PuOr   :PuRd   :Purples        :RdBu   :RdGy   :RdPu   :RdYlBu :RdYlGn :Reds   :ScalarMappable        :Set1   :Set2   :Set3   :Spectral       :Wistia :YlGn   :YlGnBu :YlOrBr :YlOrRd :afmhot        :autumn :binary :bone   :brg    :bwr    :cbook  :cividis        :cmap_d :cmaps_listed   :colors        :cool   :coolwarm       :copper :cubehelix      :datad  :flag   :get_cmap       :gist_earth   :gist_gray       :gist_heat      :gist_ncar      :gist_rainbow   :gist_stern     :gist_yarg      :gnuplot       :gnuplot2       :gray   :hot    :hsv    :inferno        :jet    :ma     :magma  :mpl    :nipy_spectral :np     :ocean  :pink   :plasma :prism  :rainbow        :register_cmap  :seismic        :spring        :summer :tab10  :tab20  :tab20b :tab20c :terrain        :turbo  :twilight       :twilight_shifted      :unregister_cmap        :viridis        :winter
#All hues can add "_r" for reverse.
#For ncolors() and grp2color(), if cmap is qualitative, N should be not larger than the available color number, or setting loop=true to loop the color order.
#See also: blendcolor, bluered
#Xiong Jieyi, March 10, 2015 > Aug 14, 2015 >Feb 3, 2016 >Dec 2, 2016>20 Dec 2019 > 1 Mar 2022>25 Aug 2022

export colormap, ncolors, grp2color
export ncolor #Disused, will be removed in future.
function ncolors(N::Integer, cmap::Union{Symbol, AbstractString}=N<10 ? :C : :rainbow; ends_gap::Bool=false, rev::Bool=false, loop::Bool=false)
    if cmap==:C
        Nmax=10
        V=if loop
            mod.(0:(N-1), Nmax)
        else
            N<=Nmax || error("For qualitative cmap $cmap, N(=$N) should <= $Nmax.")
            0:(N-1)
        end
        if rev
            reverse!(V)
        end
        matplotlib.colors.to_hex.(f"C$1".(V))
    else
        qcmap=Dict(:tab20b => 20, :Pastel2 => 8, :Set1 => 9, :tab10 => 10, :Set2 => 8, :tab20c => 20, :Paired => 12, :Pastel1 => 9, :Dark2 => 8, :Accent => 8, :tab20 => 20, :Set3 => 12)
        Nmax=get(qcmap, Symbol(cmap), nothing)
        if isnothing(Nmax)
            colormap(1:N, cmap; ends_gap=ends_gap, rev=rev)
        else
            V=if loop
                mod.(0:(N-1), Nmax)
            else
                N<=Nmax || error("For qualitative cmap $cmap, N(=$N) should <= $Nmax.")
                0:(N-1)
            end
            if rev
                reverse!(V)
            end
            matplotlib.colors.to_hex.(Base.eachrow(getproperty(PyPlot.cm, cmap)(V)))
        end
    end
end
function ncolor(N::Integer; hue::Union{Symbol, AbstractString}=N<10 ? :C : :rainbow, kw...)
    @warn("ncolor() has been renamed by ncolors() in . The old one will be abolished in future.")
    ncolors(N, hue; kw...)
end
colormap(N::Integer; hue::Union{Symbol, AbstractString}=N<10 ? :C : :rainbow)=error("colormap(N; hue=...) has been changed to ncolors(n; cmap=...) in 1 Mar 2022.")
function colormap(X::AbstractVector{<:Real}, cmap::Union{Symbol, AbstractString}=:rainbow;
                  hue::Union{Symbol, AbstractString, Nothing}=nothing, # hue=... is disused, will be removed in the future.
                  min_::Real=NaN, max_::Real=NaN,
                  rev::Bool=false, ends_gap::Bool=false)
    if !isnothing(hue)
        @warn("hue=... is abolished. Use cmap=... instead.")
        cmap=hue
    end
    if isnan(min_) || isnan(max_)
        mi, ma=extrema(X)
        gap=ends_gap ? 0.5*(ma-mi)/(length(X)-1) : 0
        if isnan(min_)
            min_=mi-gap
        end
        if isnan(max_)
            max_=ma+gap
        end
    end
    V=if rev
        max.(min.((X.-min_)./(max_-min_), 1.0), 0.0)
    else
        max.(min.((max_.-X)./(max_-min_), 1.0), 0.0)
    end
    matplotlib.colors.to_hex.(Base.eachrow(getproperty(PyPlot.cm, cmap)(V)))
end
function grp2color(G::fastgrp, cmap::Union{Symbol, AbstractString}=G.grpnum<10 ? :C : :rainbow;
                   hue::Union{Symbol, AbstractString, Nothing}=nothing,
                   kw...)
    if !isnothing(hue)
        @warn("hue=... is abolished. Use cmap=... instead.")
        cmap=hue
    end
    uC=ncolors(G.grpnum, cmap; kw...)
    C=fill("", G.itemnum)
    for i=1:G.grpnum
        l=want(G, i)
        C[l].=uC[i]
    end
    C
end
grp2color(G::Group, arg...)=grp2color(fastgrp(G), arg...)
#}}

#{{ blendcolor
#[[ blendcolor ]]
# mixed_color_hex = blendcolor(colorA, colorB[, colorA_proportion = 0.5])
# Blend two colors and output a hex string.
# See also: colormap, ncolors, grp2color
# 4 Jul 2022

export blendcolor
function blendcolor(A, B, P::Real=0.5)
    if P<0 || P>1
        error("The 3rd argument (A color ratio) should between 0 and 1.")
    end
    A=matplotlib.colors.to_rgb(A)
    B=matplotlib.colors.to_rgb(B)
    matplotlib.colors.to_hex((A.*P).+(B.*(1-P)))
end
#}}

#{{ pcaplot
#[[ pcaplot ]]
# handles = pcaplot(Matrix, ...; method=:pca|:pca_R|:tSNE_py|:mds, grp=Group|fastgrp, legend_inp=inp([labels ;]loc=:best) | legend_kw=(...)|nothing, ..., dot_inp=inp("o")|[inp("ro"), inp("bx"), ...], tsne_init_pc=50, tsne_kw=...) #Draw in dots.
# ... = pcaplot(Matrix, ...; label=..., grp=Group|fastgrp, label_colorset=[...], ...) #draw in text.
# ... = pcaplot(distance_vector|distance_matrix; method=:mds, ...)
#Draw 2D PCA, tSNE or classical MDS plot. For PCA, data will be centerred but not be scaled. For MDS, input is a distance matrix. If you only need the coordinates but not drawing, try princomp() instead.
#Each row is a sample (dot).
#legend_inp is the 2~end input arguments of PyPlot.legend function. If labels is given, the order of labels is according to the alphabetic order of 'color=Group'.
#dot_inp is the plot input or txtplot input (label != nothing). If color has been assigned, it could be a container for each color group, otherwise it should be a container for each row of data.
#label_colorset is the color set with the length of group number in grp. Must be used with label=... and grp=....
#In tSNE_py mode,
##   tsne_init_pc assign the initial PC number. By default, tsne_init_pc=50 if size(M, 2)>100, or no pre-PCA otherwise.
#    tsne_kw assigns the arguments of sklearn.manifold.TSNE. See https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html
#useR: using R:prcomp or MultivariateStats.jl to calculate PCA.
#18 Oct 2019: For history reason, key argument 'grp=' has a alias 'color='.
#See also: boxplot, density, treeplot, princomp, kdescatter
#Xiong Jieyi, 25 Sep 2014 >Feb 3, 2016>Jun 9, 2016>2 Nov 2017>3 Jun 2018>18 Oct 2019

function pcaplot(M::AbstractMatrix{TT}, arg1...; dryrun::Bool=false, color=nothing, grp=color, label=nothing, legend_inp=inp(loc=:best), legend_kw=legend_inp[2], showlegend::Bool=legend_kw!=nothing, dot_inp=inp("o"), label_colorset=nothing, method::Symbol=:pca, tsne_init_pc::Int=size(M, 2)>100 ? 50 : 0, tsne_kw=NamedTuple(), arg2...) where {TT<:Real}

    (pcaM, pc_contri)=princomp(M; method=method, tsne_init_pc=tsne_init_pc, tsne_kw=tsne_kw)
    # pp=importpkg(:PyPlot, preloaded=true)
    pp=PyPlot
    G = grp==nothing ? nothing : isa(grp,fastgrp) ? grp : fastgrp(grp)
    if label==nothing
        if G!=nothing
            (isa(dot_inp,Tuple{Tuple, Any}) || length(dot_inp)==G.grpnum) || error(f"When color has assigned, doc_inp (length=$1) should has the length equal to color unique number ($2)."(length(dot_inp)), G.grpnum)
            h=grploop(G,pcaM[:,1],pcaM[:,2],input_grpno=true) do i,xx,yy
                c_dot_inp=isa(dot_inp,Tuple{Tuple, Any}) ? dot_inp : dot_inp[i]
                pp.plot(xx,yy,c_dot_inp[1]...,arg1...;arg2...,c_dot_inp[2]...)
                
            end
            showlegend && pp.legend(G.grpid; legend_kw...)
        elseif isa(dot_inp,Tuple{Tuple, Any})
            h=pp.plot(pcaM[:,1],pcaM[:,2],dot_inp[1]...,arg1...;arg2...,dot_inp[2]...)
        else
            size(pcaM,1)==length(dot_inp) || error(f"When color has not assigned, doc_inp (length=$1) should has the length equal to data row number ($2)."(length(doc_inp)), size(pcaM,1))
            h=map(1:size(pcaM,1)) do i
                pp.plot(pcaM[i,1],pcaM[i,2],dot_inp[i][1]...,arg1...;arg2...,dot_inp[i][2]...)[1]
            end
        end
    else
        if G!=nothing
            if isa(dot_inp,Tuple{Tuple, Any}) && isempty(dot_inp[2])
                cc=isnothing(label_colorset) ? grp2color(G) : grpfunexp((i, _)->getrow(label_colorset, i), G, 1:G.itemnum, no=true)
                h=txtplot(pcaM[:,1],pcaM[:,2],label,arg1...; colors=cc, arg2...)
                rowloop(cc[G.grpfirsti,:], as_vec=true) do c
                    pp.plot(NaN,NaN,"-";color=c) #For legend
                end
            else
                h=grploop(G,pcaM[:,1],pcaM[:,2],label,input_grpno=true) do i,xx,yy,clabel
                    c_dot_inp=isa(dot_inp,Tuple{Tuple, Any}) ? dot_inp : dot_inp[i]
                    ch=txtplot(xx,yy,clabel,arg1...;arg2...,c_dot_inp[2]...)
                    pp.plot(NaN,NaN,"-";color=ch[1].get_color()) #For legend
                end
            end
            if showlegend
                if isempty(legend_inp[1])
                    pp.legend(G.grpid ; legend_kw...)
                else
                    pp.legend(legend_inp[1]...; legend_kw...)
                end
            end
        else
            h=txtplot(pcaM[:,1],pcaM[:,2],label,arg1...; arg2..., dot_inp[2]...)
        end
    end
    if pc_contri!=nothing
        pp.xlabel(f"PC1 $1%"(num2str(pc_contri[1]*100,maxdec=2)))
        pp.ylabel(f"PC2 $1%"(num2str(pc_contri[2]*100,maxdec=2)))
    end
    h
    # end
    
    #Or using Python's PCA
    #Need: @pyimport matplotlib.mlab as mlab
    #
    # opca=pycall(mlab.PCA,PyAny,M)
    # plot(Geom.point,
    #      Guide.xlabel("PC1 $(opca[:fracs][1]*100)%"),
    #      Guide.ylabel("PC2 $(opca[:fracs][2]*100)%"),
    #      arg1...;
    #      x=opca[:Y][:,1],y=opca[:Y][:,2],arg2...)
end
function pcaplot(V::AbstractVector, arg...; method::Symbol=:pca, kw...)
    method==:mds || error("When input is a (distance) vector, method should only be :mds .")
    pcaplot(distvec2mx(V), arg...; method=method, kw...)
end
export pcaplot

#}}

#{{ treeplot
#[[ treeplot ]]
#treeplot(Matrix; label=1:N, dist=Distances.Euclidean, texttree=true, column_width=80)
#Draw neighbor-join tree using BioPython.
#column_width only work when texttree is true.
#See also: pcaplot
#Xiong Jieyi, 25 Sep 2014

function needBioPhylo()
    if !isdefined(:py_Bio_Phylo_TreeConstruction)
        require("PyCall")
        @eval Main.PyCall.@pyimport Bio.Phylo as py_Bio_Phylo
        @eval Main.PyCall.@pyimport Bio.Phylo.TreeConstruction as py_Bio_Phylo_TreeConstruction
    end
end
function treeplot(M::Matrix{T};label::AbstractVector=1:size(M,1),
                  texttree::Bool=true,  column_width::Real=80, dist::DataType=Void) where {T<:Real}
    needBioPhylo()
    isdefined(:Distances)||require("Distances")
    
    if !iseltype(label,AbstractString)
        label=[string(x) for x in label]
    end
    dist=dist<:Void ? Main.Distances.Euclidean : dist
    D=Main.Distances.pairwise(Main.Distances.dist(),M')
    dismx=PyVector([PyVector(vec(D[i,1:i])) for i=1:size(M,1)])
    M=py_Bio_Phylo_TreeConstruction._DistanceMatrix(names=label,matrix=dismx)
    treeconstructor=py_Bio_Phylo_TreeConstruction.DistanceTreeConstructor()
    tree=treeconstructor.nj(M)
    tree.ladderize()
    if texttree
        py_Bio_Phylo.draw_ascii(tree, column_width=column_width)
    else
        py_Bio_Phylo.draw_graphviz(tree)
    end
end
export treeplot
#}}

#{{ kmeancostplot
#[[ kmeancostplot ]]
#kmeanscostplot(X, ks::AbstractVector; kmeans_param=(), Gadfly_param...)
#Draw the cost curve of kmeans.
#See also: grpboxplot, kmenas
#Xiong Jieyi, 2 Oct 2014

function kmeancostplot(X::Matrix, ks::AbstractVector; kmeans_param=(), wargs...)
    plot(x=ks,y=kmeancost(X,ks; kmeans_param...),Geom.point,Geom.line;wargs...)
end
export kmeancostplot
#}}

#{{ showimfile
#[[ showimfile ]]
#showimfile("filename")
#Show a bit-mapped file on screen, using PyPlot package.
#In the non-ipython model, function draw using PyPlot.
#See also: rfig, filelink, dirlink, pltsave
#Xiong Jieyi, November 9, 2014>December 6, 2014

function showimfile(filename::AbstractString,eqaxis::Bool=false)
    if !(isdefined(Main, :IJulia) && Main.IJulia.inited)
        #We can using python to draw
        isdefined(Main,:PyPlot) || require("PyPlot")
        PP=Main.PyPlot
        #PP.figure()
        eqaxis && PP.axis("equal")
        PP.xticks([],[])
        PP.yticks([],[])
        PP.box("off")
        PP.imshow(PP.imread(filename))
    else
        #Or, more simply:
        ext=match(r"(?<=\.)\w+$",filename).match
        display(MIME("image/$ext"),open(readbytes,filename))
    end
    nothing
end
export showimfile
#}}

#{{ filelink dirlink
#[[ filelink dirlink ]]
#[IJulia] filelink(filename)
#[IJulia] dirlink(path = ".")
#Show the link of a file or all the files in a given path.
#Only work in IJulia. In other environment, function only return nothing.
#See also: showimfile, rfig
#Xiong Jieyi, March 18, 2015

export filelink, dirlink
function filelink(filename::AbstractString)
    if isijulia()
        importpy(:IPython!display).FileLink(filename)
    else
        nothing
    end
end
function dirlink(filename::AbstractString=".")
    if isijulia()
        importpy(:IPython!display).FileLinks(filename)
    else
        nothing
    end
end
#}}

#{{ mousetext
#[[ mousetext ]]
#mousetext(AbstractString|Num)
#Write text to figure using mouse.
#See also:
#Xiong Jiyei, January 1, 2015

function mousetext(X)
    # plt=importpkg(:PyPlot, preloaded=true)
    plt=PyPlot
    X=string(X)
    println("Click figure for $X")
    plt.text(plt.ginput()[1]...,X;
             horizontalalignment="center", verticalalignment="center")
end
export mousetext
#}}

#{{ clustermap
#[[ clustermap ]]
#(RowTree, RowCluster, clustermapObj) = clustermap( Matrix;
#     using_fastcluster_for_row = true,
#     row_ids= c"...", col_ids= c"...", #Default: 0-based index
#     row_linkage= ... or below line:
#        row_corr= false,
#        row_method= "average"(D)|"complete",
#        row_matric= "euclidean"(D)|"correlation",
#        row_cut= cutoff|row_cutbyN= MaxClusterNumber,
#     row_grp= Group,
#     row_cluster= true, #For all above params, `row` can replaced by `col`.
#     split_color= "w",
#     label_hue= :C, #i.e. C0, C1, ...,  C9. Or use cmap like :Paired.
#     z_score|standard_scale= 0(rows)|1(columns),
#     xticklabels= T|F|vec, #Whether show x ticks or not. The default value is (Matrix_column_num<=40).
#     yticklabels= T|F|vec, #Whether show y ticks or not. The default value is (Matrix_row_num<=45).
#     figsize=(width, height), cmap=..., col_id_rotation=0  )
#Draw a clustermap using seaborn package. The row will be clustered by fastcluster python package instead to boost its speed (Can be repressed by setting using_fastcluster_for_row=false). However, the default column clustering is not changed. It is useful to first build clusters and choose cutoff, then using clustermap(..., row_linkage=xxx , cut=... ) to show visualize these clusters.
# row_corr= true only change the default row_method and row_matrix to "complete" and "correlation" respectively. Otherwise they will be "average" and "euclidean" by default.
# split_color is the color of gap to split different row groups.
#
# seaborn.clustermap: https://seaborn.pydata.org/generated/seaborn.clustermap.html
# seaborn.heatmap https://seaborn.pydata.org/generated/seaborn.heatmap.html
# seaborn.cmap see https://matplotlib.org/examples/color/colormaps_reference.html
#See also: grpheatmap, Rheatmap, colormap, fastcluster
#Xiong Jieyi, Aug 14, 2015>Aug 20, 2015>Feb 14, 2016>Feb 14, 2016>Dec 4, 2016>Jan 9, 2017>5 May 2017>2 Jan 2018>11 Dec 2020

export clustermap
function clustermap(data::AbstractMatrix{T};
                    col_ids::Vector=AbstractString[],
                    row_ids::Vector=AbstractString[],
                    using_fastcluster_for_row::Bool=true,
                    using_fastcluster_for_col::Bool=true,
                    row_cluster::Bool=true,
                    col_cluster::Bool=true,
                    row_corr::Bool=false,
                    row_method::Union{AbstractString,Symbol}=row_corr ? "complete" : "average",
                    row_metric::Union{AbstractString,Symbol}=row_corr ? "correlation" : "euclidean",
                    col_corr::Bool=false,
                    col_method::Union{AbstractString,Symbol}=col_corr ? "complete" : "average",
                    col_metric::Union{AbstractString,Symbol}=col_corr ? "correlation" : "euclidean",
                    row_cut::Real=NaN, row_cutbyN::Int=0,
                    col_cut::Real=NaN, col_cutbyN::Int=0,
                    label_hue::Union{AbstractString,Symbol}=:C,
                    row_grp::Union{Group,Void}=nothing,
                    col_grp::Union{Group,Void}=nothing,
                    row_linkage=nothing,
                    col_linkage=nothing,
                    split_color::Union{AbstractString,Symbol}="w",
                    col_id_rotation::Real=NaN,
                    xticklabels::Bool=size(data,2)<=40, #Not sure if it is a good default number.
                    yticklabels::Bool=size(data,1)<=45, #If >45, y stick labels will be stick together and seems ugly.
                    wargs...) where {T<:Real}
    
    # pp=importpkg(:PyPlot)
    pp=PyPlot
    pd=importpy(:pandas)
    # sb=importpy(:seaborn!apionly)
    sb=importpy(:seaborn)
    
    addwarg=Any[]
    if using_fastcluster_for_row
        if row_cluster && isnothing(row_linkage)
            row_linkage=fastcluster(data,method=row_method,metric=row_metric)
            push!(addwarg, (:row_linkage, row_linkage))
        end
        if !isnan(row_cut)
            row_grp=fcluster(row_linkage, row_cut)
        elseif row_cutbyN>0
            row_grp=fcluster(row_linkage, row_cutbyN, criterion="maxclust")
        end
    else
        if !(row_cluster==true && row_corr==false && string(row_method)=="average" && string(row_metric)=="euclidean" && isnan(row_cut) && row_cutbyN==0)
            error("When using_fastcluster_for_row=false, all row-clustering parameters should be in default.")
        end
        if !isnothing(row_linkage)
            push!(addwarg, (:row_linkage, row_linkage))
        end
    end
    
    if using_fastcluster_for_col
        if col_cluster && isnothing(col_linkage)
            col_linkage=fastcluster(data', method=col_method, metric=col_metric)
            push!(addwarg, (:col_linkage, col_linkage))
        end
        if !isnan(col_cut)
            col_grp=fcluster(col_linkage, col_cut)
        elseif col_cutbyN>0
            col_grp=fcluster(col_linkage, col_cutbyN, criterion="maxclust")
        end
    else
        if !(col_cluster==true && col_corr==false && string(col_method)=="average" && string(col_metric)=="euclidean" && isnan(col_cut) && col_cutbyN==0)
            error("When using_fastcluster_for_col=false, all col-clustering parameters should be in default.")
        end
        if !isnothing(col_linkage)
            push!(addwarg, (:col_linkage, col_linkage))
        end
    end
    
    if !isnothing(row_grp)
        push!(addwarg,(:row_colors,grp2color(row_grp,Symbol(label_hue))))
    end
    if !isnothing(col_grp)
        push!(addwarg,(:col_colors,grp2color(col_grp,Symbol(label_hue))))
    end

    t=Any[]
    if !isempty(row_ids)
        push!(t,(:index,row_ids))
    end
    if !isempty(col_ids)
        push!(t,(:columns,col_ids))
    end
    obj=sb.clustermap(pd.DataFrame(data;t...); row_cluster=row_cluster, col_cluster=col_cluster, yticklabels=yticklabels, xticklabels=xticklabels, addwarg..., wargs...)

    # Show right row label in horizental (The default is in vertical).
    pp.setp(obj.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    if !isnan(col_id_rotation)
        pp.setp(obj.ax_heatmap.xaxis.get_majorticklabels(), rotation=col_id_rotation)
    end

    if row_grp!=nothing && !isempty(row_grp)
        if row_cluster
            L=leaves(row_linkage)
        else
            L=1:rownum(row_grp)
        end
        ax=obj.ax_heatmap
        for i=1:length(L)-1
            if getrow(row_grp,L[i])!=getrow(row_grp,L[i+1])
                # ax[:axhline](length(L)-i, c=split_color) # Used in seaborn <0.8.1
                ax.axhline(i, c=split_color) # Used since seaborn 0.8.1
            end
        end
    end

    (row_grp,row_linkage,obj)
end
#}}

#{{ grpheatmap
#[[ grpheatmap ]]
#clustermapObj = grpheatmap( Matrix; row_grp=Group, row_cluster=true,
#     row_ids=c"...", col_ids=c"...",
#     row_corr=false,
#     row_method="average"|"complete", row_matric="euclidean"|"correlation", cut= cutoff|cutbyN= MaxClusterNumber,
#     col_grp=Group, label_hue="Paired", split_color="w"
#     z_score|standard_scale=0(rows)|1(columns), col_cluster=true, figsize=..., )
#Draw cluster map for each group respectively. Data will first be devided by row_grp, then do the cluster in each group, and draw heatmap in a one figure.
# row_corr=true only change the default row_method and row_matrix to "complete" and "correlation" respectively. Otherwise they will be "average" and "euclidean" by default.
#If row_grp or col_grp are segregated, the label and element number of each group will be shown.
# Problem solve: In the case some fonts cannot be showed in pdf file by pltsave, try importpy("seaborn").set(font="Bitstream Vera Sans") before drawing.
#See also: clustermap
#Xiong Jieyi, Aug 24, 2015 > Jan 9, 2017 > 5 May 2017

export grpheatmap
function grpheatmap(data::AbstractMatrix{T};
                    row_ids::Vector=AbstractString[],
                    row_grp::Union{Group, Nothing}=nothing,
                    row_corr::Bool=false,
                    row_method::Union{AbstractString,Symbol}=row_corr ? "complete" : "average",
                    row_metric::Union{AbstractString,Symbol}=row_corr ? "correlation" : "euclidean",
                    row_cluster::Bool=true,
                    col_grp::Union{Group, Nothing}=nothing,
                    wargs...) where {T<:Real}
    if !isnothing(row_grp)
        idx=Int[]
        grploop(fastgrp(row_grp, stable=true), data, collect(1:size(data,1))) do D, I
            if row_cluster && length(I)>1
                Z=fastcluster(D,method=row_method,metric=row_metric)
                append!(idx,I[leaves(Z)])
            else
                append!(idx,I)
            end
        end
        data=data[idx,:]
        row_grp=getrow(row_grp,idx)
        if !isempty(row_ids)
            row_ids=getrow(row_ids,idx)
        end
        row_cluster=false
    end

    o=clustermap(data;row_grp=row_grp,row_ids=row_ids,row_cluster=row_cluster,col_grp=col_grp,wargs...)[3]
    
    #Add row group label
    if !isnothing(row_grp) && issegregatedr(row_grp)
        # (txty,num),unigrp=grpfun(x->(mean(x),length(x)),row_grp,rownum(row_grp)-collect(1:rownum(row_grp))+1,multi_output=true) #Used in seaborn <0.8.1
        (txty,num),unigrp=grpfun(x->(mean(x),length(x)),row_grp,collect(1:rownum(row_grp)),multi_output=true) #Used since seaborn 0.8.1
        for i=1:length(txty)
            o.ax_row_colors.text(0,txty[i],string(getrow(unigrp,i))*"\nN=$(num[i])",rotation="vertical",ha="center",va="center")
        end
    end

    #Add column group label
    if !isnothing(col_grp) && issegregatedr(col_grp)
        (txtx,num),unigrp=grpfun(x->(mean(x),length(x)),col_grp,collect(1:rownum(col_grp)),multi_output=true)
        for i=1:length(txtx)
            o.ax_col_colors.text(txtx[i],1,string(getrow(unigrp,i))*"\nN=$(num[i])",ha="center",va="center")
        end
    end
        
    o
end

#}}

#{{ venn
#[[ venn ]]
# venn(V1, V2[, V3]; label=c"A,B,C")
# Draw 2-sets or 3-sets venn diagram using matplotlib_venn package.
# V1...3 could be ID vectors or Bool vectors. For Bool vectors, all inputs should be in the same lengths.
# Tip: Set the font size:
# o=venn(...)
# for text in [o.set_labels; o.subset_labels]
#     text.set_fontsize(14)
# end
# # .subset_labels are the numbers on the circles.
# See also: setcmpr
#Xiong Jieyi, Jan 14, 2016 > Sep 19, 2016

export venn
function venn(A::AbstractVector{T}, B::AbstractVector{T}; label=c"A,B") where {T}
    # importpkg(:PyPlot)
    ob=importpy(:matplotlib_venn)
    Base.invokelatest(ob.venn2, (pyfun(:set,A),pyfun(:set,B)), label)
end
function venn(A::AbstractVector{T}, B::AbstractVector{T}, C::AbstractVector{T}; label=c"A,B,C") where {T}
    # importpkg(:PyPlot)
    ob=importpy(:matplotlib_venn)
    ob.venn3((pyfun(:set,A),pyfun(:set,B),pyfun(:set,C)),label)
end
function venn(Xs::AbstractVector{Bool}...;wargs...)
    N=length(Xs[1])
    all(x->length(x)==N,Xs[2:end]) || error("Bool-vector input should be in the same lengths.")
    L=(1:N)
    venn(map(x->L[x],Xs)...;wargs...)
end
#}}

#{{ PDFfigs
#[[ PDFfigs next finish ]]
#pdf=PDFfigs(filename; subplot=(RowNum, ColNum), figsize=(W, H), tight_layout=false, savefig_kw=(...,))
#for i=1:N
#   row_i, column_i = next(pdf; newpage=false) #Note that running "next" is required before drawing the first figure.
#   #Draw the i-th figure here ...
#end
#finish(pdf; dpi=...)
#Draw multipage PDF figures. Support subplot. By default, figure will be saved as bbox_inches="tight".
#If tight_layout is true, each figure will be PyPlot.tight_layout() before saving.
#If newpage=true, the remaining empty subplots in this page will be given up.
#If dpi is given, the output figure will be rasterized by external software in the last step (Need 'convert' be installed).
#savefig_kw is the parameters of matplotlib.backends.backend_pdf.savefig(). e.g., savefig_kw=(dpi=300,) will change the DPI when rasterized=true in the axes.
#See also: pltsave, filelink, dirlink
#Xiong Jieyi, March 24, 2015 >Dec 23, 2016>14 Jan 2019>29 Jul 2019>9 Dec 2019

export PDFfigs, next, finish
mutable struct PDFfigs
    pdf
    subplot::Tuple{Int,Int}
    curplot::Int
    figsize::Union{Void,Tuple{Int,Int}}
    tight_layout::Bool
    filename::AbstractString
    savefig_kw
    
    function PDFfigs(filename::AbstractString; subplot::Tuple{Int,Int}=(1,1), figsize::Union{Void,Tuple{Int,Int}}=nothing, tight_layout::Bool=false, savefig_kw=NamedTuple())
        if !ismatch(r"\.pdf$",filename)
            filename=filename*".pdf"
        end
        # plt=importpkg(:PyPlot)
        plt=PyPlot
        bknd=importpy(:matplotlib!backends!backend_pdf)
        pdf=Base.invokelatest(bknd.PdfPages, filename)
        new(pdf,subplot,0,figsize,tight_layout,filename,savefig_kw)
    end
end
function next(o::PDFfigs,newpage::Bool=false)
    # plt=importpkg(:PyPlot)
    plt=PyPlot
    if newpage || o.curplot==0 || o.curplot>=prod(o.subplot)
        if o.curplot>0
            o.tight_layout &&  plt.tight_layout()
            o.pdf.savefig(; bbox_inches="tight", o.savefig_kw...)
            plt.close()
            o.curplot=0
        end
        o.figsize==nothing || plt.figure(figsize=o.figsize)
    end
    o.curplot+=1
    if o.subplot!=(1,1)
        plt.subplot(o.subplot...,o.curplot)
    end
    (div(o.curplot-1, o.subplot[2])+1, rem(o.curplot, o.subplot[2]))
end
function finish(o::PDFfigs; dpi::Int=-1)
    o.tight_layout && tight_layout()
    o.pdf.savefig(; bbox_inches="tight", o.savefig_kw...)
    o.pdf.close()
    if dpi>0
        println("Rasterizing vector graphic to $dpi dpi...")
        figfile=o.filename
        run(`convert -density $dpi -compress Zip $figfile $figfile`)
    end
end
#}}

#{{ pdf_start pdf_next pdf_finish
#[[ pdf_start pdf_next pdf_finish ]]
# `Using PDFfigs type is more recommended than these functions.`
# pdf_start(pdf_filename = "~/.tempfig.pdf")
# PyPlot.pygui(false) #Optional.
# #Draw page 1 using PyPlot here...
# pdf_next(is_close_fig=true) #Next-page
# #Draw page 2 using PyPlot here...
# pdf_next() #Next-page
# ...
# pdf_next() #To save the last page, pdf_next() should be ran again before pdf_finish(.)
# pdf_finish() #Finish the whold pdf files.
#
#Draw a multiple-page pdf figure.
#See also: PDFfigs, pltsave, rfig, showimfile, filelink, dirlink
#Xiong Jieyi, March 18, 2015 >Feb 3, 2016

export pdf_start, pdf_next, pdf_finish
function pdf_start(figfile::AbstractString="/home/jxiong/.tempfig.pdf")
    global PYPLOT_MULTIPAGE_PDF
    @assert(!isdefined(:PYPLOT_MULTIPAGE_PDF) || isa(PYPLOT_MULTIPAGE_PDF,Void),
            "Last pdf file is unfinished. Try pdf_finish() first.")
    # importpkg(:PyPlot) #A error will occur if PyPlot is not added in advance.
    bknd=importpy(:matplotlib!backends!backend_pdf)
    if !ismatch(r"\.pdf$",figfile)
        figfile=figfile*".pdf"
    end
    PYPLOT_MULTIPAGE_PDF=bknd.PdfPages(figfile)
    nothing
end
function pdf_next(is_close_fig::Bool=true)
    global PYPLOT_MULTIPAGE_PDF
    @assert(isdefined(my,:PYPLOT_MULTIPAGE_PDF) && !isa(PYPLOT_MULTIPAGE_PDF,Void),
            "Have not ran pdf_start(...) yet.")
    # bknd=importpy(:matplotlib!backends!backend_pdf)
    PYPLOT_MULTIPAGE_PDF.savefig()
    if is_close_fig
        # plt=importpkg(:PyPlot)
        plt=PyPlot
        plt.close()
    end
    nothing
end
function pdf_finish()
    global PYPLOT_MULTIPAGE_PDF
    @assert(isdefined(my,:PYPLOT_MULTIPAGE_PDF) && !isa(PYPLOT_MULTIPAGE_PDF,Void),
            "Have not ran pdf_start(...) yet.")
    PYPLOT_MULTIPAGE_PDF.close()
    PYPLOT_MULTIPAGE_PDF=nothing
    # plt=importpkg(:PyPlot)
    # plt.pygui(true)
    nothing
end

#}}

#{{ guideline
#[[ guideline ]]
# handle_vec=guideline("-|/"[, intercept=0.0|[...]], linesty="r:";
# label="labels of each line"|labels=StringVector, text_kw=(horizontalalignment="right", verticalalignment="bottom", fontsize=10), label_loc=1|2(default).
# PyPlot.plot_params...)
#Draw a guideline. Default is draw "-|/".
#intercept only work in -(y) or |(x) line, but not for / and anti-slash line. For convenience, it can also be a vector for parallel lines.
#If label is given, a label will marked on the left/bottom (label_loc=1) or right/top (label_loc=2) end of line.
#"labels=Vector" only works when intercept is also a vector.
#See also: fitline, freelegend, tight_layout
#Xiong Jieyi, March 25, 2015 > Aug 28, 2016 > Jan 9, 2017

export guideline
function guideline(typ::AbstractString="-|/",B::Real=0.0,linesty::AbstractString="r:";
                   label::String="", text_inp::Tuple{Tuple, Any}=inp(), text_kw=text_inp[2], label_loc::Integer=2, wargs...)
    # plt=importpkg(:PyPlot)
    plt=PyPlot
    xlim=plt.xlim()
    ylim=plt.ylim()
    h=plt.PyObject[]
    if '-' in typ
        append!(h,plt.plot(xlim,[B,B],linesty;wargs...))
        isempty(label) || plt.text(xlim[label_loc], B, label; horizontalalignment=["left","right"][label_loc], verticalalignment="bottom", fontsize=10, text_kw...)
    end
    if '|' in typ
        append!(h,plt.plot([B,B],ylim,linesty;wargs...))
        isempty(label) || plt.text(B, ylim[label_loc], label; horizontalalignment="right", verticalalignment=["bottom","top"][label_loc], rotation="vertical", fontsize=10, text_kw...)
    end
    if '/' in typ
        lim=[max(xlim[1],ylim[1]),min(xlim[2],ylim[2])]
        append!(h,plt.plot(lim,lim,linesty;wargs...))
        isempty(label) || plt.text(lim[label_loc], lim[label_loc], label; horizontalalignment=["left","right"][label_loc], verticalalignment=["bottom","top"][label_loc], rotation=45, fontsize=10, text_kw...)
    end
    if '\\' in typ
        lim=[max(xlim[1],ylim[1]),min(xlim[2],ylim[2])]
        append!(h,plt.plot([lim[2],lim[1]],lim,linesty;wargs...))
        isempty(label) || plt.text(lim[label_loc], reverse(lim, dims=1)[label_loc], label; horizontalalignment=["left","right"][label_loc], verticalalignment=["top","bottom"][label_loc], rotation=-45, fontsize=10, text_kw...)
    end
    plt.xlim(xlim)
    plt.ylim(ylim)
    vcat(h...)
end
guideline(typ::AbstractString, linesty::AbstractString; wargs...)=guideline(typ,0.0,linesty;wargs...)
function guideline(typ::AbstractString, Bs::AbstractVector{T}, args...; labels::Vector=String[], wargs...) where {T<:Real}
    if isempty(labels)
        vcat(map(B->guideline(typ,B,args...;wargs...),Bs)...)
    else
        vcat(map((B,L)->guideline(typ,B,args...;label=L, wargs...),Bs, labels)...)
    end
end
#}}

#{{ fitline
#[[ fitline ]]
# h, m, b = fitline(X, Y[, "-."]; across0=false, xends=(min, max), plot_params...)
# Draw a 1D fit line. The xlim and ylim will not be changed by fitline().
# across0: whether the line across point (0, 0);
# To make the line across the whole figure, set xends=xlim().
#See also: guildline, cortest
#Xiong Jieyi, 27 Oct 2017 > 2 May 2018

export fitline
function fitline(X::AbstractVector,Y::AbstractVector,linsty::AbstractString="-.",args...; deg::Int=1, xends=nothing, across0::Bool=false, wargs...)
    # pp=importpkg(:PyPlot, preloaded=true)
    pp=PyPlot
    l=.!isnan.(X) .& .!isnan.(Y) .& .!isinf.(X) .& .!isinf.(Y)
    m,b=if across0
        sum(l)>=1 || error("No valid value.")
        (pp.plt.np.linalg.lstsq(reshape(X[l], :, 1), Y[l])[1], 0.0)
    else
        sum(l)>=2 || error("Should be at least two valid values.")
        pp.plt.np.polyfit(X[l], Y[l], 1)
    end
    if xends==nothing
        xends=[minimum(X[l]), maximum(X[l])]
    elseif isa(xends, Tuple)
        xends=vcat(xends...)
    end
    xlm=pp.xlim()
    ylm=pp.ylim()
    h=pp.plot(xends, m.*xends+b, linsty, args...; wargs...)
    pp.xlim(xlm)
    pp.ylim(ylm)
    (h, m, b)
end
fitline(X::Real, Y::Real, args...; wargs...)=fitline([zero(X), X], [zero(Y), Y], args...; wargs...)
#}}

#{{ kdescatter
#[[ kdescatter ]]
# (DensityValue, handle) = kdescatter (X, Y, ...; s=6, cmap="jet", densityfun=..., normto=(min, max), scatter_kw...)
# Or DensityValue = kdescatter(X, Y; draw=false) #No draw.
# Draw scatter plot with the color of density. The parameter is identical to scatter function, except don't set c=.... The kde output will be processed by the given densityfun.
# normto=(m1, m2) will re-scale density value to interval [m1, m2], and set scatter(...; vmin=0, vmax=1, ...) if draw=true. e.g.: cmap="binary", normto=(0.2, 1) will draw dots from gray to black.
# See also: pcaplot, grpplot
#Xiong Jieyi, 18 Apr 2018 > 12 Nov 2021

export kdescatter
function kdescatter(X::AbstractVector, Y::AbstractVector; draw::Bool=true, densityfun::Function=x->x, s=6, cmap="jet", normto::Union{Tuple, Nothing}=nothing, wargs...)
    l=.!(isnan.(X) .| isnan.(Y) .| isinf.(X) .| isinf.(Y))
    if any(.!l)
        @warn f"$1 data are either NaN or Inf. They are ignored."(sum(.!l))
    end
    @pkgfun(KernelDensity, kde, InterpKDE, pdf)
    ob=InterpKDE(kde((X[l], Y[l])))
    C=fill(NaN, length(X))
    C[l]=densityfun(map((x, y)->pdf(ob, x, y), X[l], Y[l]))
    if !isnothing(normto)
        Cmin, Cmax=extrema(C[l])
        C[l]=((C[l].-Cmin)./(Cmax-Cmin)).*(normto[2]-normto[1]).+normto[1]
    end
    if draw
        @pkgfun(PyPlot, scatter)
        h=if isnothing(normto)
            scatter(X, Y; c=C, s=s, cmap=cmap, wargs...)
        else
            scatter(X, Y; c=C, s=s, cmap=cmap, vmin=0, vmax=1, wargs...)
        end
        (C, h)
    else
        C
    end
end
#}}

#{{ retangular
#[[ retangular ]]
# h = retangular((x1,y1), (x2,y2), ...)
#Draw a retangular on the figure.
#Parameters:
#  fill (default: true); color; edgecolor or ec; facecolor or fc; label; linestyle or ls (like: "none" ":", "-."); linewidth or lw; zorder
#See also: zebra
#Xiong Jieyi, 2 May 2018

export retangular
function retangular(P1::Tuple, P2::Tuple; wargs...)
    # plt=importpkg(:PyPlot, preloaded=true)
    plt=PyPlot
    P0=(min(P1[1], P2[1]), min(P1[2], P2[2]))
    W=abs(P1[1]-P2[1])
    H=abs(P1[2]-P2[2])
    rect=plt.plt.matplotlib.patches.Rectangle(P0, W, H; wargs...)
    plt.gca().add_patch(rect)
end
#}}

#{{ zebra
#[[ zebra ]]
# zebra(boundaries; hori=false, onfirst=false, other_retangular_param...)
# Draw gray strips (facecolor="0.1") on the background (zorder=0).
# hori: horizental strips (true) or vertical strips (false).
# onfirst: starts from strip or empty.
# See also: retangular
# Xiong Jieyi, 5 Mar 2019

export zebra
function zebra(X; hori::Bool=false, onfirst::Bool=false, kw...)
    @pkgfun(PyPlot, xlim, ylim)
    xm=xlim()
    ym=ylim()
    offset=onfirst ? 0 : 1
    if hori
        X=[ym[1]; X[ym[1].<X.<ym[2]]; ym[2]]
        for i=(1+offset):2:(length(X)-1)
            retangular((xm[1], X[i]), (xm[2], X[i+1]); zorder=0, fill=true, ec="none", fc="0.9", kw...)
        end
    else
        X=[xm[1]; X[xm[1].<X.<xm[2]]; xm[2]]
        for i=(1+offset):2:(length(X)-1)
            retangular((X[i], ym[1]), (X[i+1], ym[2]); zorder=0, fill=true, ec="none", fc="0.9", kw...)
        end
    end
end
#}}

#{{ tight_suptitle
#[[ tight_suptitle ]]
# tight_suptitle("title"; topgap=0.05, suptitle_param...)
# tight_suptitle(;topgap=0.05)  #Keep top empty place while tight_layout others.
# Draw super title and apply the tight_layout(). topgap (0-1) is the relative space for the super title.
# See also: freelegend, guideline
# Xiong Jieyi, 18 Jun 2018

export tight_suptitle
function tight_suptitle(arg...; topgap::Real=0.05, warg...)
    # plt=importpkg(:PyPlot, preloaded=true)
    plt=PyPlot
    h=isempty(arg) ? nothing : plt.suptitle(arg...; warg...)
    plt.tight_layout(rect=[0, 0, 1, 1-topgap])
    h
end
#}}

#{{ freelegend
#[[ freelegend ]]
# freelegend(c"r.-, 0.5x, b[], ..."|"ro", c"label1, label2, label3, ..."; each_kw=(mfc=["r", missing, ...], ...), all_kw=...; legend_args...)
# Or: freelegend([inp(..), inp(...), ...], c"label1, label2, ..."; legend_args...)
# freelegend()  # Remove current legend.
# Draw legend freely or remove it.
# The each_kw is by element (x[i]), or by row as a vector (x[i, :]) only when input is a matrix.
# each_kw element can also be missing, which will ignore this keyword argument.
# First input like "0.5o" is allowed in freelengend() although it is not allowed in plot(). This function just call plot(..., "o", color="0.5", ...).
# Input can be "[]" or "r[]", where function will  use matplotlib.patches.Rectangle() instead plot() to draw, and the color code before ("r") will be assigned to facecolor (color="r").
# Some useful legend parameters:
# loc=
# # 'best'		0
# # 'upper right'	1
# # 'upper left'	2
# # 'lower left'	3
# # 'lower right'	4
# # 'right'		5
# # 'center left'	6
# # 'center right'	7
# # 'lower center'	8
# # 'upper center'	9
# # 'center'		10
# Tip: Set loc=6, bbox_to_anchor=(1.02, 0.5) will make the legend at right outside.
# ncol=1 : The number of columns that the legend has.
# other legend parameters: title, title_fontsize, fancybox=false.
# For detail see https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.legend.html
# See also: tight_layout, guideline, kdescatter, retangular
# Xiong Jieyi, 18 Jun 2018 > 25 Aug 2018 > 4 Jul 2021

export freelegend
function freelegend(marks, labels::AbstractVector; each_inp::Tuple{Tuple, Any}=inp(), each_kw=each_inp[2], all_inp::Tuple{Tuple, Any}=inp(), all_kw=all_inp[2],  wargs...)
    # plt=importpkg(:PyPlot, preloaded=true)
    plt=PyPlot
    h=map(1:length(labels)) do i
        m=isa(marks, AbstractVector) ? marks[i] : marks
        kwargs=collect(Pair{Symbol, Any}, pairs(all_kw))
        if isa(m, Tuple{Tuple, Any})
            append!(kwargs, m[2])
            m=m[1][1]
        end
        for (k, v) in pairs(each_kw)
            cv=if isa(v, AbstractMatrix)
                v[i, :]
            else
                v[i]
            end
            if !ismissing(cv)
                push!(kwargs, k => cv)
            end
        end
        if endswith(m, "[]")
            t=if length(m)>2
                kwargs=[:color=>m[1:end-2], kwargs...]
            end
            plt.gca().add_patch(plt.matplotlib.patches.Rectangle((NaN, NaN), NaN, NaN; kwargs...))
        else
            if !isnothing((t=match(r"^(0\.\d+)(.*)$", m);))
                kwargs=[:color=>t[1], kwargs...]
                m=t[2]
            end
            plt.plot(NaN, NaN, m; kwargs...)[1]
        end
    end
    plt.legend(h, labels; wargs...)
end
function freelegend()
    # plt=importpkg(:PyPlot, preloaded=true)
    plt=PyPlot
    h=plt.gca().get_legend()
    if h!=nothing
        h.remove()
    end
end
#}}

#{{ ax_height ax_width ax_rl_x ax_rl_y
#[[ ax_height ax_width ax_rl_x ax_rl_y ]]
# (p *) axis_height = ax_height(p=1.0)
# (p *) axis_width  = ax_width(p=1.0)
# relative_X = ax_rl_x(p=0.5)
# relative_Y = ax_rl_y(p=0.5)
#Quickly access axis height, width and relative coordinates.
#18 Feb 2022

export ax_height, ax_width, ax_rl_x, ax_rl_y
function ax_height(p::Real=1)
    # plt=importpkg(:PyPlot, preloaded=true)
    plt=PyPlot
    a, b=plt.ylim()
    (b-a)*p
end
function ax_width(p::Real=1)
    # plt=importpkg(:PyPlot, preloaded=true)
    plt=PyPlot
    a, b=plt.xlim()
    (b-a)*p
end
function ax_rl_x(p::Real=0.5)
    # plt=importpkg(:PyPlot, preloaded=true)
    plt=PyPlot
    a, b=plt.xlim()
    (b-a)*p+a
end
function ax_rl_y(p::Real=0.5)
    # plt=importpkg(:PyPlot, preloaded=true)
    plt=PyPlot
    a, b=plt.ylim()
    (b-a)*p+a
end
#}}

addhelpfromfile(@__FILE__, inmodule=@__MODULE__)
end
