@info "loading and compiling packages..."


using Revise, Plots, Statistics, Plots.PlotMeasures, JLD, PrettyTables, Dates, CSV
using DataFrames
using DSP


# custom pkgs:
# julia_path = joinpath(ENV["HOME"],"dev/nanoseis/julia_src")
# push!(LOAD_PATH,julia_path)
# #push!(LOAD_PATH,joinpath(julia_path,"HDF5utilities"))
# push!(LOAD_PATH,joinpath(julia_path,"HDF5utilities/src"))
# push!(LOAD_PATH,joinpath(julia_path,"PumpCurveTools/src"))

using PumpCurveTools

using HDF5utilities

using DMD


ENV["JULIA_DEBUG"]=""

#printstyled("debug environmental variable $(ENV["JULIA_DEBUG"])=\n",color=:magenta)

##  2 


if !isdir(fig_dir)
    dirlist = splitpath(fig_dir)
    if length(dirlist) > 1
        for idx_dir in 1:length(dirlist)-1
            if !isdir(dirlist[idx_dir])
                mkdir(dirlist[idx_dir])
            end
        end
    end
    mkdir(fig_dir)
end


if !@isdefined do_clusters
    do_clusters=true
end

if !@isdefined company
    company = "halliburton"
end

# if !@isdefined expt_date_list
#     expt_date_list = map(x->date_from_name(x,company=company), infileroot_list)
# end


# if !@isdefined pumpcurve_file
#     pumpcurve_df = DataFrame()
# elseif pumpcurve_file == "none"
#     pumpcurve_df = DataFrame()
# else
#     pumpcurve_df = load_pumpcurve_data(pumpcurve_file,time_shift=pumpcurve_time_shift, company=company)
#     @info "***** 1st pump data time=$(pumpcurve_df[1,:truedatetime])"
# end

if do_production
    production_df = CSV.File(production_file) |> DataFrame!
end

if do_clusters
    cluster_df = CSV.File(cluster_file) |> DataFrame!;
end
function parse_comma(s::String)
    parts = split(s,",")
    return parse(Float64,parts[1])*1000 + parse(Float64,parts[2])
end
"""
    get_cluster_mds_hal(df::DataFrame, idx_row::Int)

For a Halliburton file, get cluster MDs for the given row `idx_row`, which  
corresponds to stage.  For that stage, extract MDs for each cluster in that stage. 
Returns `Array{Float64,1}` with those values.
"""
function get_cluster_mds_hal(df::DataFrame, idx_row::Int)
    cluster_mds = [parse_comma(df[idx_row,3])]
    for i in 4:18
        if typeof(df[idx_row,i]) == Missing
            break
        end
        push!(cluster_mds,parse_comma(df[idx_row,i]))
    end
    return cluster_mds
end
"""
    get_cluster_mds_schlum(df::DataFrame, idx_cluster::Int)

For cluster `idx_cluster`, extract reference MD for each cluster (average of
provided top and bottom MD).
"""
function get_cluster_mds_schlum(df::DataFrame, idx_cluster::Int)
    # we have either 4 or 5 clusters per stage
    cluster_depths = Array{Float64,2}(undef,5,2)
    n_clusters = 0
    for i in 1:5
        top_td = df[idx_cluster, "Cluster $(i) Top TD"]
        if top_td != 0
            bot_td = df[idx_cluster, "Cluster $(i) Bottom TD"]
            cluster_depths[i,1] = top_td
            cluster_depths[i,2] = bot_td
            n_clusters += 1
        end
    end
    return vec(mean(cluster_depths,dims=2))
end


stage_labels, measured_signals, md_values, tm_values, tm_abs_values = DMD.read_measures_mds_times(input_dir, idx_level);

if do_clusters
    cluster_md_list = Array{Array{Float64,1},1}(undef,length(stage_labels))
    for i in 1:length(stage_labels)
        if company == "halliburton"
            cluster_md_list[i] = get_cluster_mds_hal(cluster_df,parse(Int,stage_labels[i]))
        elseif company == "schlumberger"
            cluster_md_list[i] = get_cluster_mds_schlum(cluster_df,parse(Int,stage_labels[i]))
        end
    end
end

# note: tm_values has times relative to original data set start, and the values match those displayed in extracted mode heatmaps

# make plots before proppant time window
plots_combined = Plots.plot();
for idx_stage in 1:length(get_stages)

    # idx_stage = 1

    # find index in data read from disk by read_measures_mds_times, as that data will have all stages computed
    local idx_stage_raw = findfirst(x->x==get_stages[idx_stage],stage_labels)

    # if findfirst(x->x==stage_labels[idx_stage], get_stages) != nothing
    local idx_pre_proppant = findfirst(x->x>=pre_proppant_times[idx_stage], tm_values[idx_stage_raw])

    if (abs(tm_values[idx_stage_raw][idx_pre_proppant-1]-pre_proppant_times[idx_stage])
        < abs(tm_values[idx_stage_raw][idx_pre_proppant]-pre_proppant_times[idx_stage]) )
        idx_pre_proppant -= 1
    end

    local mds_unstaggered = [(md_values[idx_stage_raw][i]+md_values[idx_stage_raw][i+1])/2 for i in 1:length(md_values[idx_stage_raw])-1]
    mds_unstaggered .*= md_scalar
    if zero_ref_mds
        local md_shift = mds_unstaggered[1]
        mds_unstaggered .-= md_shift
    end

    if do_clusters
        cluster_mds = copy(cluster_md_list[idx_stage_raw])
        cluster_mds .-= md_shifts[idx_stage]
        cluster_mds .*= md_scalar
        if zero_ref_mds
            cluster_mds .-= md_shift
        end
    end

    if zero_ref_mds
        local dist_label = "Distance ($(md_units))"
    else
        local dist_label = "Measured depth ($(md_units))"
    end

    local p = plot(mds_unstaggered, measured_signals[idx_stage_raw][:,idx_pre_proppant],label="",color=:gray,
             xlabel=dist_label,ylabel="Acoustic energy indicator",
                   xlimits = (xlimits_list[idx_stage][1],xlimits_list[idx_stage][2]),right_margin=3mm, left_margin=3mm)
    global plots_combined
    plot!(plots_combined,mds_unstaggered, measured_signals[idx_stage_raw][:,idx_pre_proppant],label="$(stage_labels[idx_stage])",
         xlabel=dist_label,ylabel="Acoustic energy indicator",
         xlimits = (xlimits_list[idx_stage][1],xlimits_list[idx_stage][2]),right_margin=3mm, left_margin=3mm)
    
    if do_clusters
        y_cluster = fill(cluster_plot_y_pre[idx_stage],length(cluster_mds))
        scatter!(p, cluster_mds, y_cluster, markersize=3, label="",markerstrokewidth=0, color=:black)
        scatter!(plots_combined, cluster_mds, y_cluster, markersize=3, label="",markerstrokewidth=0, color=:black)
    end
    if do_production
        plot!(p, production_df[:,:DepthM], production_df[:,"Zonal Gas Rate"]/500.0, markersize=3, label="",markerstrokewidth=0, color=:black)
    end

    @info "pre proppant code: writing file $(joinpath(fig_dir,"mode_cluster_level$(idx_level)_stage$(get_stages[idx_stage])_preprop.png"))"
    savefig(p,joinpath(fig_dir,"mode_cluster_level$(idx_level)_stage$(get_stages[idx_stage])_preprop.png"))

end
if do_production
    plot!(plots_combined, production_df[:,:DepthM], production_df[:,"Zonal Gas Rate"]/500.0, markersize=3, label="",markerstrokewidth=0, color=:black)
end
savefig(plots_combined, joinpath(fig_dir,"mode_clusters_all_level$(idx_level)_preprop.png"))

# make plots after proppant time window
plots_combined = Plots.plot();
for idx_stage in 1:length(get_stages)

    # idx_stage = 2

    # find index in data read from disk by read_measures_mds_times, as that data will have all stages computed
    local idx_stage_raw = findfirst(x->x==get_stages[idx_stage],stage_labels)

    # if findfirst(x->x==stage_labels[idx_stage], get_stages) != nothing
    idx_post_proppant = findfirst(x->x>=post_proppant_times[idx_stage], tm_values[idx_stage_raw])

    if (abs(tm_values[idx_stage_raw][idx_post_proppant-1]-post_proppant_times[idx_stage])
        < abs(tm_values[idx_stage_raw][idx_post_proppant]-post_proppant_times[idx_stage]) )
        idx_post_proppant -= 1
    end

    local mds_unstaggered = [(md_values[idx_stage_raw][i]+md_values[idx_stage_raw][i+1])/2 for i in 1:length(md_values[idx_stage_raw])-1]
    mds_unstaggered .*= md_scalar

    if do_clusters
        cluster_mds = copy(cluster_md_list[idx_stage_raw])
        cluster_mds .-= md_shifts[idx_stage]
        cluster_mds .*= md_scalar
    end

    if zero_ref_mds
        local md_shift = mds_unstaggered[1]
        mds_unstaggered .-= md_shift
        if do_clusters
            cluster_mds .-= md_shift
        end
    end

    if zero_ref_mds
        local dist_label = "Distance ($(md_units))"
    else
        dist_label = "Measured depth ($(md_units))"
    end

    local p = plot(mds_unstaggered, measured_signals[idx_stage_raw][:,idx_post_proppant],label="",color=:gray,
             xlabel=dist_label,ylabel="Acoustic energy indicator",
             xlimits = (xlimits_list[idx_stage][1],xlimits_list[idx_stage][2]),right_margin=3mm, left_margin=3mm)
    global plots_combined
    plot!(plots_combined,mds_unstaggered, measured_signals[idx_stage_raw][:,idx_post_proppant],label="$(stage_labels[idx_stage])",
         xlabel=dist_label,ylabel="Acoustic energy indicator",
         xlimits = (xlimits_list[idx_stage][1],xlimits_list[idx_stage][2]),right_margin=3mm, left_margin=3mm)
    
    if do_clusters
        y_cluster = fill(cluster_plot_y_post[idx_stage],length(cluster_mds))
        scatter!(p, cluster_mds, y_cluster, markersize=3, label="",markerstrokewidth=0, color=:black)
        scatter!(plots_combined, cluster_mds, y_cluster, markersize=3, label="",markerstrokewidth=0, color=:black)
    end
    if do_production
        plot!(p, production_df[:,:DepthM], production_df[:,"Zonal Gas Rate"]/500.0, markersize=3, label="",markerstrokewidth=0, color=:black)
    end

    savefig(p,joinpath(fig_dir,"mode_cluster_level$(idx_level)_stage$(get_stages[idx_stage])_postprop.png"))

end
if do_production
    plot!(plots_combined, production_df[:,:DepthM], production_df[:,"Zonal Gas Rate"]/500.0, markersize=3, label="",markerstrokewidth=0, color=:black)
end
savefig(plots_combined, joinpath(fig_dir,"mode_clusters_all_level$(idx_level)_postprop.png"))
