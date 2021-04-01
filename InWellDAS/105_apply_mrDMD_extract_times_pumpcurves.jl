#=

apply mrDMD to time windows in a single data sets

key difference - in this approach, we extract data time range from HDF5 files, and get pump start time
from pump curve file.  Given a user-specified dτ, we then compute t0_dmd and t1_dmd to get as much of
data file as possible

Paramters must specify time window length dτ, data file name, pump curve file etc

=#

@info "executing apply_mrDMD"

## 0.5 load code

@info "loading and compiling packages..."

using Revise, Plots, Statistics, Plots.PlotMeasures, JLD, Dates, CSV, DataFrames
# using DSP

# custom pkgs:
julia_path = joinpath(ENV["HOME"],"dev/nanoseis/julia_src")
push!(LOAD_PATH,julia_path)
push!(LOAD_PATH,joinpath(julia_path,"SignalUtilities"))
push!(LOAD_PATH,joinpath(julia_path,"HDF5utilities/src"))
push!(LOAD_PATH,joinpath(julia_path,"PumpCurveTools/src"))
#using DASFVReader, SignalUtilities

using PumpCurveTools

using HDF5utilities

include(joinpath(julia_path,"SignalUtilities/src/SignalUtilities.jl"))
include(joinpath(julia_path,"DASFVReader/src/DASFVReader.jl"))
include(joinpath(julia_path,"DMD/src/DMD.jl"))
#include(joinpath(julia_path,"HDF5utilities/HDF5utilities.jl"))


#include("../9999_primary_flows/201_pumpcurve_tools.jl")



# use first line to turn on debugging printouts, will get messy and slow though...
#ENV["JULIA_DEBUG"]=Main
ENV["JULIA_DEBUG"]=""



## 1 - load data  #################################################

@info "loading data and applying mrDMD......"

#cd(@__DIR__)

if !isdir(fig_dir)
    mkdir(fig_dir)
end
if !isdir(out_dir)
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
datetimestamp,tm_vals_list = HDF5utilities.get_datetimestamp_timevalues(data_filename,time_dim=1);
date_end = Date(datetimestamp)
if tm_vals_list[end] < tm_vals_list[1]
    date_end += Day(1)
end
# bounding DateTimes for the dataset
datetimeBnds = (datetimestamp, DateTime(date_end, Time(Hour(tm_vals_list[end]),
                                                       Minute(tm_vals_list[end]),
                                                       Second(tm_vals_list[end])
                                                       )
                                        )
                )
data_length_seconds = (datetimeBnds[2]-datetimeBnds[1]).value/1000 # seconds

@time pumpcurve_df = load_pumpcurve_data(pumpcurve_file,time_shift=time_shift)

@info "Pump curve data summary:"
println(describe(pumpcurve_df))
println("\n\n")

pumpcurve_df_forstage = filter(row->(row.truedatetime>=datetimeBnds[1] && row.truedatetime <=datetimeBnds[2]), pumpcurve_df);
function abs2reltime(dt::DateTime, refdt::DateTime)
    t0 = Time(refdt).instant.value/10^9
    t1 = Time(dt).instant.value/10^9
    if Date(dt) != Date(refdt)
        t1 += 24*60*60
    end
    t1 - t0
end
pumpcurve_df_forstage = hcat(pumpcurve_df_forstage, DataFrame(:reldatetime=>map(x->abs2reltime(x,pumpcurve_df_forstage[1,:truedatetime]), pumpcurve_df_forstage[!,:truedatetime])))
@info "Pump curve data summary for current stage only:"
println(describe(pumpcurve_df_forstage))
println("\n\n")

pump_start_idx = get_pump_start_index(pumpcurve_df_forstage, treat_pressure_threshold)
@info "pumping starts at pump curve data frame index $(pump_start_idx)"
# get number of seconds of data before pumping starts
prepump_time_seconds = datetime2unix(pumpcurve_df_forstage[pump_start_idx,:truedatetime]) - datetime2unix(datetimeBnds[1])
prepump_time_seconds = floor(Int,prepump_time_seconds/dτ_dmd)*dτ_dmd
lead_time_seconds = min(max_time_before_pump_seconds,prepump_time_seconds)
@info "Check prepump_time_seconds=$(prepump_time_seconds), lead_time_seconds=$(lead_time_seconds)"
t0_dmd = round(Int,(prepump_time_seconds - lead_time_seconds)/dt)*dt
if t0_dmd<dt  # make sure we don't try access sample at index 0
    t0_dmd = dt
end
t_post_t0 = data_length_seconds - t0_dmd
n_dmd = floor(Int,t_post_t0/dτ_dmd) - 1
t1_dmd = t0_dmd + n_dmd*dτ_dmd

@info "Bounding DateTimes of data are: $(datetimeBnds[1]), $(datetimeBnds[2])\n\tTrying to access data $(lead_time_seconds) seconds before pump onset"
@info "Set initial mrDMD application time to $(t0_dmd) seconds, final time is $(t1_dmd)"


# these values bound the true reference locations:
d_md = md_values[2] - md_values[1] # assume constant
md_values_staggered = [md_values[i] - d_md/2 for i in 1:length(md_values)]
md_values_staggered = vcat(md_values_staggered,[md_values[end] + d_md/2])
md_values_staggered .*= 3.28

i0_dmd = round(Int,t0_dmd/dt)
t0_dmd_abs = tm_vals_list[i0_dmd]
i1_dmd = round(Int,t1_dmd/dt)
n_per_dmd = round(Int,dτ_dmd/dt)

i0 = 1
i1 = size(data,1) # last time sample index
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
p_datadmd = heatmap(ts_abs_dmd,md_values,data[i0_dmd:di:i1_dmd,:]',fill=(true,cmap_data),clims=(-clip,clip),ylabel="Measured depth (feet)",xlabel="Time ($(time_units))",yflip=true,right_margin=0mm)
savefig(joinpath(fig_dir,"data_dmdrange.png"))




## 2 - apply multiwindow dmd ##################################################
# this does a transpose as well as reduce size
data_input = transpose(data[i0_dmd:i1_dmd,:]);
# plot data subset for sanity check, compare to pump - can reuse tm_abs_dmd from above
ts_rel_dmd = map(x->( (x.instant.value)/10^9 - ts_abs_dmd[1].instant.value/10^9), ts_abs_dmd)
p_datadmd = heatmap(ts_rel_dmd,md_values .* 3.28,data_input[:,1:di:end],fill=(true,cmap_data),clims=(-clip,clip),ylabel="Measured depth (feet)",xlabel="Time ($(time_units))",yflip=true,right_margin=5mm,colorbar=:none,xlimits=(0,t1_dmd-t0_dmd));
p_pump = plot(pumpcurve_df_forstage[:,:reldatetime],pumpcurve_df_forstage[:,"Treating Pressure"],label="",ylabel="P",right_margin=5mm,xlimits=(0,t1_dmd-t0_dmd));
l = @layout [a{0.2h} ; b]
p = plot(p_pump, p_datadmd, layout=l)

# @time dmd_results = apply_dmd2array(data_input, ,i0_dmd, i1_dmd, ndmd_modes)
@time dmd,dmdgrid_sc = DMD.apply_mrdmd2windows(data_input,Float32(dτ_dmd),Float32(dt),ndmd_modes,ndmd_levels,subtract_modes=false, include_dc=true);

p_datadmd = heatmap(ts_rel_dmd .+ t0_dmd, md_values .- md_values[1], data_input[:,1:di:end], fill=(true,cmap_data), clims=(-clip,clip),
                    ylabel="Distance (m)",xlabel="Time ($(time_units))",yflip=true,right_margin=5mm,colorbar=:none,
                    xlimits=(t0_dmd,t1_dmd),dpi=150);
p_dmd = DMD.heatmap_mrdmd_grid(dmd,dmdgrid_sc,Float32(t0_dmd),Float32(dt),title="",cmap=cmap_dmd)
min_dmd,max_dmd = extrema(dmdgrid_sc)
plot(p_datadmd, plot(p_dmd,clims=(0,max_dmd/5),colorbar=:none,right_margin=5mm),layout=(2,1),dpi=150);
savefig(joinpath(fig_dir,"dmd_$(dτ_dmd)sec.png"))

p_datadmd = heatmap(ts_rel_dmd .+ t0_dmd, md_values .- md_values[1], data_input[:,1:di:end], fill=(true,cmap_data), clims=(-clip,clip),
                    ylabel="Distance (m)",xlabel="Time ($(time_units))",yflip=true,right_margin=5mm,colorbar=:none,
                    xlimits=(t0_dmd,t1_dmd),dpi=300);
p_dmd = DMD.heatmap_mrdmd_grid(dmd,dmdgrid_sc,Float32(t0_dmd),Float32(dt),cmap=:grays)
min_dmd,max_dmd = extrema(dmdgrid_sc)
plot(p_datadmd, plot(p_dmd,clims=(0,max_dmd/5),colorbar=:none,right_margin=5mm),layout=(2,1),dpi=300);
savefig(joinpath(fig_dir,"dmd_$(dτ_dmd)sec_anon.png"))



## 3 - extract max amplitude mode display #################################

@info "generating images of max modes..."

# plot max amplitudes modes for all levels
pm_list = Array{Plots.Plot,1}(undef,0)
for i in 1:length(dmd[1])
    # i = 1
    global clip, pm_list
    ts,chans,maxmodegrid = DMD.mrDMD_modegrid(dmd,i,min_chan_idx=channel_list[1],t0=Float64(t0_dmd))
    ts_abs = [t0_dmd_abs + Dates.Millisecond(round(Int,(ts[i]-ts[1])*1000.0)) for i in 1:length(ts)]
    #heatmap(ts,chans,maxmodegrid,clims=(-2,2),fill=(true,cmap_data),xlabel="Time (sec)",ylabel="Channel")
    modestd = std(maxmodegrid)
    clip=modestd*2
    pm = heatmap(ts_abs,md_values_staggered,maxmodegrid,clims=(-clip,clip),fill=(true,cmap_dmd),xlabel="Time (UTC)",ylabel="Measured depth (feet)",yflip=true)
    push!(pm_list,pm)
    savefig(pm,joinpath(fig_dir,"maxmodegrid_$(i).png"))
end


## 4 - Apply deglitch to generate final mode estimates - also take absolute values #################

@info "deglitching modes to get final estimates..."

median_filter_length = 11
median_threshold=2.5
pmd_list = Array{Plots.Plot,1}(undef,length(dmd[1]))
for i in 1:length(dmd[1])
     # i = 1
    @info "working on level $(i)...."
    @info "\tt0_dmd=$(t0_dmd)"
    global clip
    ts,chans,maxmodegrid = DMD.mrDMD_modegrid(dmd,i,min_chan_idx=channel_list[1],t0=Float64(t0_dmd))
    maxmodegrid_deglitch = SignalUtilities.deglitch_data(abs.(maxmodegrid),median_filter_length,threshold_factor=median_threshold);
    # set up absolute UTC time values
    ts_abs = [t0_dmd_abs + Dates.Millisecond(round(Int,(ts[i]-ts[1])*1000.0)) for i in 1:length(ts)]
    pmd_list[i] = heatmap(ts_abs, md_values_staggered, (maxmodegrid_deglitch),clims=(0,3),title="Mode: level $(i)",xlabel="Time (UTC)", ylabel="Measured depth (feet)",yflip=true,colorbar=:none);
    local l = @layout [a{0.2h} ; b]
    local p = plot(p_pump, pmd_list[i], layout=l)
    savefig(p,joinpath(fig_dir,"mode_level$(i)_deglitched.png"))
    @info "\tsaving mrDMD results..."
    jldopen(joinpath(out_dir,outfileroot*"_maxmode_deglitch_$(i).jld"), "w") do file
        write(file, "modegrid_deglitch", maxmodegrid_deglitch)
    end
    @info "\t\tCheck min and max time: $(ts[1])   $(ts[end])"
    jldopen(joinpath(out_dir,outfileroot*"_maxmode_times_$(i).jld"), "w") do file
        write(file, "modegrid_times", ts)
    end
    jldopen(joinpath(out_dir,outfileroot*"_maxmode_times_abs_$(i).jld"), "w") do file
        write(file, "modegrid_times_abs", ts_abs) 
    end
    jldopen(joinpath(out_dir,outfileroot*"_maxmode_chans_$(i).jld"), "w") do file
        write(file, "modegrid_chans", chans)
    end
    jldopen(joinpath(out_dir,outfileroot*"_maxmode_md_values_stagger_$(i).jld"),"w") do file
        write(file, "modegrid_md_values_stagger", md_values_staggered)
    end
end
plot(pmd_list...,layout=(3,2),size=(1000,700));
savefig(joinpath(fig_dir,"all_deglitched_modes.png"))
