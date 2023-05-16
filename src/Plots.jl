
@userplot ScatterGrid

"""
    $(SIGNATURES)

Add documentation.
"""
function scattergrid end


# pt = collect(node(g))
set_forplot(G) = eachrow(reduce(hcat, nodes(G)))

@recipe function f(sg::ScatterGrid;
                   r = 2,
                   R = 20,
                   f_ms = i -> r + R*(.5)^i, # Marker size
                   f_ma = (d, dmin, dmax) -> (d - dmin)/(dmax - dmin), # Mark alpha/opacity
                   color = "rgb(238,37,35)"
                   )
    G, = sg.args
    Dict_depth = grid_groupbydepth(G)
    I = sort!(collect(keys(Dict_depth))) # Sort it for the legend
    dmin = first(I)
    dmax = last(I)
    # Set default values
    axis --> ([-1.1,1.1],)
    # xaxis --> ("x",)
    # yaxis --> ("y",)
    legend --> :outerright
    legendfontsize --> 10
    thickness_scaling --> 1
    shift = dmin -1
    # Plot each depth group
    for d in I
        x, y = set_forplot(Dict_depth[d]) # We only plot in two dimensions.
        @series begin
            seriestype := :scatter
            label := "Depth $d"
            # markershape := :circle
            ms := f_ms(d - shift) # Marker size
            mc := color # Marker color
            ma := f_ma(d, dmin, dmax) # Mark alpha/opacity
            msw := .3 + (1 - f_ma(d, dmin, dmax))/2 # Border width (0 = No border)
            msc := color # Border color
            msa := f_ma(d, dmin, dmax) # Border alpha
            x, y
        end
    end
end
