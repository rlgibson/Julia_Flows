
using Dates, Plots, DataFrames, CSV

function load_pumpcurve_data(filename::String; time_shift=Hour(0))

    if !isfile(filename)
        error("load_pumpcurve_data: requested file $(filename) is not a file...")
    end

    df = CSV.File(filename) |> DataFrame!

    insertcols!(df, :truedatetime => map(x->DateTime(x,"m/d/y H:M:S.s") + time_shift, df[:, :textdatetime]))
end

# function pumpcurve_for_timerange(df::DataFrame, column_hdr::Symbol, tmin::DateTime, tmax::DateTime, n_timestep::Int)
#     filter(row->(row.truedatetime<=tmax && row.truedatetime>=tmin), df)[1:n_timestep:end,column_hdr]
# end


## YIKES - make sure this is accessible in DMD_analysis <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

"""
    pumpcurve_for_timerange(df::DataFrame, column_hdr::String, tmin::DateTime, tmax::DateTime, n_timestep::Int)

Return two lists, one of time samples and one of requested data values
"""
function pumpcurve_for_timerange(df::DataFrame, column_hdr::String, tmin::DateTime, tmax::DateTime, n_timestep::Int)

    tms = filter(row->(row.truedatetime<=tmax && row.truedatetime>=tmin), df)[1:n_timestep:end,:truedatetime]
    curve = filter(row->(row.truedatetime<=tmax && row.truedatetime>=tmin), df)[1:n_timestep:end,column_hdr]

    return tms,curve
end


function get_pump_start_index(pumpcurve_df::DataFrame, pressure_threshold::Float64)
    idx_stg_start = nothing

    for i in 1:size(pumpcurve_df,1)
        if pumpcurve_df[i,"Treating Pressure"] >= pressure_threshold
            #            @info "Found stage pumping start: i=$(i), t=$(pumpcurve_df[i,:truedatetime]), P=$(pumpcurve_df[i,"Treating Pressure"])"
            idx_stg_start = i
            break
        end
    end
    
    return idx_stg_start
end
