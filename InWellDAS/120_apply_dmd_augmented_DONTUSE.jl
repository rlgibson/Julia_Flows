@info "executing apply_mrDMD"

## 0.5 load code

@info "loading and compiling packages..."

using Revise, Plots, Statistics, Plots.PlotMeasures, JLD, Dates
# using DSP

# custom pkgs:

using PumpCurveTools
using SignalUtilities
using HDF5utilities
using DMD



# use first line to turn on debugging printouts, will get messy and slow though...
#ENV["JULIA_DEBUG"]=Main
ENV["JULIA_DEBUG"]=""



## 1 - load data  #################################################

@info "loading data and applying mrDMD......"

#cd(@__DIR__)

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
if !isdir(out_dir)
    dirlist = splitpath(out_dir)
    if length(dirlist) > 1
        for idx_dir in 1:length(dirlist)-1
            if !isdir(dirlist[idx_dir])
                mkdir(dirlist[idx_dir])
            end
        end
    end
    mkdir(out_dir)
end

# get raw data 
data_filename = joinpath(datadir,infile)
hdf5_id,dataset = HDF5utilities.get_hdf5_dataset(data_filename);
dt = HDF5utilities.get_rowinterval(dataset)
dx = HDF5utilities.get_colinterval(dataset)
time_units = HDF5utilities.get_rowunits(dataset)
channel_units = HDF5utilities.get_colunits(dataset)
data = HDF5utilities.get_dataarray(dataset);
channel_list = Int.(HDF5utilities.get_channel_indices(dataset));
timestamp,tm_vals_list = HDF5utilities.get_timestamp_timevalues(dataset,1);
md_values =  HDF5utilities.get_md_values(dataset);
expt_date  = HDF5utilities.get_date(dataset)
HDF5utilities.close_hdf5(hdf5_id)


# these values bound the true reference locations:
d_md = md_values[2] - md_values[1] # assume constant
md_values_staggered = [md_values[i] - d_md/2 for i in 1:length(md_values)]
md_values_staggered = vcat(md_values_staggered,[md_values[end] + d_md/2])
md_values_staggered .*= 3.28

i0_dmd =round(Int,t0_dmd/dt)
t0_dmd_abs = tm_vals_list[i0_dmd]
i1_dmd = round(Int,t1_dmd/dt)
n_per_dmd = round(Int,dτ_dmd/dt)
@info "starting index for dmd: $(i0_dmd), end index=$(i1_dmd)"

i0 = 1
i1 = size(data,1) # last time sample index
@info "\tdata array has $(size(data,1)) time samples, $(size(data,2)) channel samples"
di = 200
time_vals = [i*dt for i in 0:size(data,1)];
#channel_vals = [Int(idx_chan0 + i*dx) for i in 0:size(data,2)-1];
datastd = std(data)
clip = datastd/2
p_data = heatmap(time_vals[i0:di:i1],channel_list,data[i0:di:i1,:]',fill=(true,cmap_data),clims=(-clip,clip),ylabel="Channel",xlabel="Time ($(time_units))",yflip=true);
savefig(joinpath(fig_dir,"data_overview.png"))

di=50
# set up absolute UTC time values
ts_abs_dmd = [timestamp + Dates.Millisecond(round(Int,(time_vals[i]-time_vals[1])*1000.0)) for i in i0_dmd:di:i1_dmd]
p_datadmd = heatmap(ts_abs_dmd,md_values,data[i0_dmd:di:i1_dmd,:]',fill=(true,cmap_data),clims=(-clip,clip),ylabel="Measured depth (feet)",xlabel="Time ($(time_units))",yflip=true);
savefig(joinpath(fig_dir,"data_dmdrange.png"))



## 2 - apply multiwindow dmd ##################################################
#data_input = transpose(data[i0_dmd:i1_dmd,:]);
#data_input = transpose(data);
data_input = transpose(data[i0_dmd:i1_dmd,:]);
@time dmd,dmdgrid_sc = DMD.apply_mrdmd2windows(data_input,Float32(dτ_dmd),Float32(dt),ndmd_modes,ndmd_levels,subtract_modes=false, include_dc=true);

n_modes = 0
for idx_level in 1:length(dmd[1])
    global n_modes = Array{Int,1}(undef,length(dmd)*2^(idx_level-1))
    @info "setting number of modes list to length $(length(n_modes))"
    idx_nm = 1
    for idx_dmd in 1:length(dmd)
        for idx_pane in 1:2^(idx_level-1)
            n_modes[idx_nm] = size(dmd[idx_dmd][idx_level].pane_results[idx_pane].dmd.Φ,2)
            idx_nm += 1
        end
    end
    scatter(n_modes,markersize=2,xlabel="DMD time step index", ylabel="Mode count",label="Level $(idx_level)",legend=:topleft)
    savefig(joinpath(fig_dir,"mode_count_level$(idx_level).png"))
end


growth_vals = Array{ Array{ eltype(data),1} }( undef, length(dmd[1]))
for idx_level in 1:length(dmd[1])
    #    global growth_vals = Array{eltype(data),1}(undef,length(dmd)*2^(idx_level-1))
    global growth_vals
    gv_level = Array{eltype(data),1}(undef,length(dmd)*2^(idx_level-1))
    idx_gv = 1
    for idx_dmd in 1:length(dmd)
        for idx_pane in 1:2^(idx_level-1)
            if length(dmd[idx_dmd][idx_level].pane_results[idx_pane].dmd.λ) > 0
                gv_level[idx_gv] = real(dmd[idx_dmd][idx_level].pane_results[idx_pane].dmd.λ[end])
            else
                gv_level[idx_gv] = 0.0
            end
            idx_gv += 1
        end
    end
    growth_vals[idx_level] = gv_level
end
tmvals = [dτ_dmd/(2^(length(dmd[1])-1))*i for i in 1:length(growth_vals[length(dmd[1])])]
growth_vals[length(dmd[1])]
p = scatter(tmvals, growth_vals[length(dmd[1])], label="Level $(length(dmd[1]))",markerstrokewidth=0.001,markersize=2,legend=:topleft,xlabel="Time (sec)", ylabel="Growth rate")
#p = scatter([dτ_dmd*i for i in 1:length(growth_vals[1])],growth_vals[1],xlabel="DMD time step index", ylabel="Growth rate",label="Level $(1)",legend=:bottomright,markerstrokewidth=0.001,markersize=2)
for idx_level in length(dmd[1])-1:-1:1
    if idx_level==1
        msize = 3
        scatter!(p,[dτ_dmd/(2^(idx_level-1))*i for i in 1:length(growth_vals[idx_level])],growth_vals[idx_level],label="Level $(idx_level)",markerstrokewidth=0.001,markersize=msize,color=:red)
    else
        msize = 1
        scatter!(p,[dτ_dmd/(2^(idx_level-1))*i for i in 1:length(growth_vals[idx_level])],growth_vals[idx_level],label="Level $(idx_level)",markerstrokewidth=0.001,markersize=msize)
    end
end
p
savefig(p,joinpath(fig_dir,"growth_rate_all.png"))
# n_print = 0
# for idx_applydmd in i0_dmd:n_per_dmd:i1_dmd
#     global n_print
#     if n_print%10 == 0
#         n_print = 0
#     end
#     n_print += 1
#     dmd_results = DMD.DMD_augmented(data_input,ndmd_modes,ndmd_stacks,i0_dmd,i0_dmd+n_per_dmd-1,Float32(dt));
    
# end


# # color scale for heatmap set below...
# #p_dmd = DMD.heatmap_mrdmd_grid(dmd,dmdgrid_sc,t0_dmd,Float32(dt),title="dτ=$(dτ_dmd)",cmap=cmap_dmd);
# p_dmd = DMD.heatmap_mrdmd_grid(dmd,dmdgrid_sc,t0_dmd,Float32(dt),cmap=cmap_dmd);

# min_dmd,max_dmd = extrema(dmdgrid_sc)
# # here is where we can set color scale for mrDMD heatmap...
# if !@isdefined mrdmd_overview_clim_scale
#     mrdmd_overview_clim_scale = 5 # historic value that "worked"
# end
# plot(p_datadmd, plot(p_dmd,clims=(0,max_dmd/mrdmd_overview_clim_scale)),layout=(2,1),dpi=150, right_margin=5mm, left_margin=5mm,
#      bottom_margin=5mm, top_margin=5mm);
# savefig(joinpath(fig_dir,"dmd_$(dτ_dmd)sec.png"))



# ## 3 - extract max amplitude mode display #################################

# @info "generating images of max modes..."

# # plot max amplitudes modes for all levels
# pm_list = Array{Plots.Plot,1}(undef,0)
# for i in 1:length(dmd[1])
#     # i = 1
#     global clip, pm_list
#     ts,chans,maxmodegrid = DMD.mrDMD_modegrid(dmd,i,min_chan_idx=channel_list[1],t0=Float64(t0_dmd))
#     ts_abs = [t0_dmd_abs + Dates.Millisecond(round(Int,(ts[i]-ts[1])*1000.0)) for i in 1:length(ts)]
#     #heatmap(ts,chans,maxmodegrid,clims=(-2,2),fill=(true,cmap_data),xlabel="Time (sec)",ylabel="Channel")
#     modestd = std(maxmodegrid)
#     clip=modestd*2
#     pm = heatmap(ts_abs,md_values_staggered,maxmodegrid,clims=(-clip,clip),fill=(true,cmap_dmd),xlabel="Time (sec)",ylabel="Measured depth (feet)",yflip=true)
#     push!(pm_list,pm)
#     savefig(pm,joinpath(fig_dir,"maxmodegrid_$(i).png"))
# end


# ## 4 - Apply deglitch to generate final mode estimates - also take absolute values #################

# @info "deglitching modes to get final estimates..."

# median_filter_length = 11
# median_threshold=2.5
# pmd_list = Array{Plots.Plot,1}(undef,length(dmd[1]))
# for i in 1:length(dmd[1])
#     # i = 1
#     @info "working on level $(i)...."
#     global clip
#     ts,chans,maxmodegrid = DMD.mrDMD_modegrid(dmd,i,min_chan_idx=channel_list[1],t0=Float64(t0_dmd))
#     maxmodegrid_deglitch = SignalUtilities.deglitch_data(abs.(maxmodegrid),median_filter_length,threshold_factor=median_threshold);
#     # set up absolute UTC time values
#     ts_abs = [t0_dmd_abs + Dates.Millisecond(round(Int,(ts[i]-ts[1])*1000.0)) for i in 1:length(ts)]
#     # @info "level $(i), start time=$(ts_abs[1]), end time=$(ts_abs[end])"
#     if !@isdefined deglitched_mode_clim_scale
#         clim_deglitch_top = 3  # historic old hardwired value
#     else
#         clim_deglitch_top = maximum(maxmodegrid_deglitch)/deglitched_mode_clim_scale
#     end
#     pmd_list[i] = heatmap(ts_abs, md_values_staggered, (maxmodegrid_deglitch),clims=(0,clim_deglitch_top),title="Mode: level $(i)",xlabel="Time (sec)", ylabel="Measured depth (feet)",yflip=true);
#     savefig(pmd_list[i],joinpath(fig_dir,"mode_level$(i)_deglitched.png"))
#     @info "\tsaving mrDMD results..."
#     jldopen(joinpath(out_dir,outfileroot*"_maxmode_deglitch_$(i).jld"), "w") do file
#         write(file, "modegrid_deglitch", maxmodegrid_deglitch)
#     end
#     @info "\t\tCheck min and max time: $(ts[1])   $(ts[end])"
#     jldopen(joinpath(out_dir,outfileroot*"_maxmode_times_$(i).jld"), "w") do file
#         write(file, "modegrid_times", ts)
#     end
#     jldopen(joinpath(out_dir,outfileroot*"_maxmode_times_abs_$(i).jld"), "w") do file
#         write(file, "modegrid_times_abs", ts_abs) 
#     end
#     jldopen(joinpath(out_dir,outfileroot*"_maxmode_chans_$(i).jld"), "w") do file
#         write(file, "modegrid_chans", chans)
#     end
#     jldopen(joinpath(out_dir,outfileroot*"_maxmode_md_values_stagger_$(i).jld"),"w") do file
#         write(file, "modegrid_md_values_stagger", md_values_staggered)
#     end
# end


# # plot(pmd_list...,layout=(3,2),size=(1000,700));
# # savefig(joinpath(fig_dir,"all_deglitched_modes.png"))
