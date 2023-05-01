using PyPlot

function button_press_event(event)
    dragging[] = true
    click_coord .= [event[:xdata],event[:ydata]]
end

function button_release_event(event)
    dragging[] = false
end

function motion_notify_event(event)
    if dragging[]
        #@show "move", [event[:xdata],event[:ydata]]
        new_coord = [event[:xdata],event[:ydata]]
        xy .= new_coord .- click_coord .+ xy
        click_coord .= new_coord
        set_mask!(config,mask,xy,(u,v))
    end
end

function fluidplot_setup(config,mask,p,(u,v))
    sz = size(mask)
    cx = config.Δxy[1] * (0:sz[1]-1)
    cy = config.Δxy[2] * (0:sz[2]-1)

    dpi = 96
    fig = figure(figsize=(900/dpi,300/dpi),dpi=dpi)
    fig.canvas.mpl_connect("button_press_event",button_press_event)
    fig.canvas.mpl_connect("button_release_event",button_release_event)
    fig.canvas.mpl_connect("motion_notify_event",motion_notify_event)

    axis("equal")
    readline()
    return fig
end

function fluidplot(fig,cx,cy,mask,p,(u,v),n)
    unstagger_vel((u,v),(u_r,v_r))
    clf()
    p_ = copy(p)
    p_[mask .== 0] .= NaN
    u_r[mask .== 0] .= NaN
    v_r[mask .== 0] .= NaN

    pcolormesh(cx,cy,p_')
    #clim(-50,50); colorbar()
    ri = 1:5:size(u_r,1)
    rj = 1:5:size(u_r,2)

    quiver(cx[ri],cy[rj],u_r[ri,rj]',v_r[ri,rj]',scale=100)
    title("Timestep $n")
    axis("equal")
    sleep(0.01)
end
