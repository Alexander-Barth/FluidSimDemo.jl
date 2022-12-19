using GLMakie, GeometryBasics


#const fluidplot_setup = let

    function plotting(config,mask,p::AbstractArray{T,N},(u,v); plot_every = 1) where {T,N}
        dragging = Ref(false)
        click_coord = zeros(2)
        xy = [0.,0.]

        sz = size(mask)
        cx = config.h * (0:sz[1]-1)
        cy = config.h * (0:sz[2]-1)

        u_r = zeros(T,sz)
        v_r = zeros(T,sz)
        xy .= config.xy

        fig = Figure(resolution=(800,400))
        title = Observable("")

        ax = Axis(fig[1,1], xrectzoom=false, yrectzoom=false,
                  title = title,
                  aspect=DataAspect())

        unstagger_vel((u,v),(u_r,v_r))
        p_ = copy(p)
        p_[mask .== 0] .= NaN
        u_r[mask .== 0] .= NaN
        v_r[mask .== 0] .= NaN

        pressure = Observable(p)
        f = heatmap!(ax,cx,cy,pressure)
        ri = 1:5:size(u_r,1)
        rj = 1:5:size(u_r,2)

        scale = 0.01

        a_x = Observable([cx[i] for i in ri, j in rj][:])
        a_y = Observable([cy[j] for i in ri, j in rj][:])
        a_u = Observable([scale*u_r[i,j] for i in ri, j in rj][:])
        a_v = Observable([scale*v_r[i,j] for i in ri, j in rj][:])

        arrows!(ax,a_x,a_y,a_u,a_v)

        deregister_interaction!(ax, :rectanglezoom)

        register_interaction!(ax, :my_interaction2) do event::MouseEvent, axis
            if event.type === MouseEventTypes.leftdown
                dragging[] = true
                click_coord .= event.data
            elseif event.type === MouseEventTypes.leftup
                if dragging[]
                    new_coord = event.data
                    xy .= new_coord .- click_coord .+ xy
                    click_coord .= new_coord
                    set_mask!(config,mask,xy,(u,v))
                end
                dragging[] = false
            elseif event.type === MouseEventTypes.leftdrag
                if dragging[]
                    new_coord = event.data
                    xy .= new_coord .- click_coord .+ xy
                    click_coord .= new_coord
                    set_mask!(config,mask,xy,(u,v))
                end
            end
        end

        display(fig)
        #readline()

        function myplot(mask,p,(u,v),n)
            if (n % plot_every != 0)
                return
            end
            sleep(0.01)
            unstagger_vel((u,v),(u_r,v_r))
            p_ = copy(p)
            p_[mask .== 0] .= NaN
            u_r[mask .== 0] .= NaN
            v_r[mask .== 0] .= NaN
            pressure[] = p_
            title[] = "Timestep $n"
            #f = heatmap!(ax,cx,cy,p_)
            #f = heatmap!(ax,p_)
            #@show extrema(p)
            #clim(-50,50); colorbar()
            ri = 1:5:size(u_r,1)
            rj = 1:5:size(u_r,2)

            scale = 0.01
            a_x[] = [cx[i] for i in ri, j in rj][:]
            a_y[] = [cy[j] for i in ri, j in rj][:]
            a_u[] = [scale*u_r[i,j] for i in ri, j in rj][:]
            a_v[] = [scale*v_r[i,j] for i in ri, j in rj][:]

            display(fig)
            return nothing
        end

        return myplot
    end
#end
