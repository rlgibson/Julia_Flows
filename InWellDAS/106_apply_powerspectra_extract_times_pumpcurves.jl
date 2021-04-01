#=

compute power spectra for time windows in a single data sets - with low or bandpass filter applied 

key difference - in this approach, we extract data time range from HDF5 files, and get pump start time
from pump curve file.  Given a user-specified dτ, we then compute t0_spectra and t1_spectra to get as much of
data file as possible

Paramters must specify time window length dτ, data file name, pump curve file etc

=#



## 0 - load code ################################################################

@info "loading and compiling packages..."

using Revise, Plots, Statistics, Plots.PlotMeasures, JLD, Dates, CSV, DataFrames
using FFTW
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



## 1 - load data and pump curves  #################################################

@info "loading data......"

#cd(@__DIR__)

if !isdir(fig_dir)
    mkdir(fig_dir)
end
if !isdir(out_dir)
    mkdir(out_dir)
end


# get raw data 
@info "read data..."

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


tm_vals_seconds_rel = vcat([0.0], [(tm_vals_list[i]-tm_vals_list[1]).value/10^9 for i in 2:length(tm_vals_list) ])

tm_vals_seconds_rel_staggered = vcat( tm_vals_seconds_rel .- dt/2, [tm_vals_seconds_rel[end]+dt/2])



datetimestamp = DateTime(expt_date, tm_vals_list[1])
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

# pump curve values just for the current stage

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

# now identify time range of interest from treating pressure values,
# and set time bounding values for signal analysis

pump_start_idx = get_pump_start_index(pumpcurve_df_forstage, treat_pressure_threshold)
@info "pumping starts at pump curve data frame index $(pump_start_idx)"
# get number of seconds of data before pumping starts
prepump_time_seconds = datetime2unix(pumpcurve_df_forstage[pump_start_idx,:truedatetime]) - datetime2unix(datetimeBnds[1])
prepump_time_seconds = floor(Int,prepump_time_seconds/dτ_dmd)*dτ_dmd
lead_time_seconds = min(max_time_before_pump_seconds,prepump_time_seconds)
@info "Check prepump_time_seconds=$(prepump_time_seconds), lead_time_seconds=$(lead_time_seconds)"
t0_spectra = round(Int,(prepump_time_seconds - lead_time_seconds)/dt)*dt
if t0_spectra<dt  # make sure we don't try access sample at index 0
    t0_spectra = dt
end
t_post_t0 = data_length_seconds - t0_spectra
n_spectra = floor(Int,t_post_t0/dτ_window) - 1
t1_spectra = t0_spectra + n_spectra*dτ_window

@info "Bounding DateTimes of data are: $(datetimeBnds[1]), $(datetimeBnds[2])\n\tTrying to access data $(lead_time_seconds) seconds before pump onset"
@info "Set initial spectrum time to $(t0_spectra) seconds, final time is $(t1_spectra)"


# these values bound the true reference locations:
d_md = md_values[2] - md_values[1] # assume constant
md_values_staggered = [md_values[i] - d_md/2 for i in 1:length(md_values)]
md_values_staggered = vcat(md_values_staggered,[md_values[end] + d_md/2])
md_values_staggered .*= 3.28

i0_spectra = round(Int,t0_spectra/dt)
t0_spectra_abs = tm_vals_list[i0_spectra]
i1_spectra = round(Int,t1_spectra/dt)
n_per_spectra = round(Int,dτ_window/dt)

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
ts_abs_spectra = [Time(datetimestamp) + Dates.Millisecond(round(Int,(time_vals[i]-time_vals[1])*1000.0)) for i in i0_spectra:di:i1_spectra]
p_data4spectra = heatmap(ts_abs_spectra,md_values,data[i0_spectra:di:i1_spectra,:]',fill=(true,cmap_data),clims=(-clip,clip),ylabel="Measured depth (feet)",xlabel="Time",yflip=true,right_margin=0mm)
savefig(joinpath(fig_dir,"data_spectrarange.png"))

p_data4spectra_rel = heatmap(time_vals[i0_spectra:di:i1_spectra],md_values./3.28084 .- md_values[1]/3.28084,data[i0_spectra:di:i1_spectra,:]',fill=(true,cmap_data),clims=(-clip,clip),ylabel="Distance (m)",xlabel="Time ($(time_units))",yflip=true,right_margin=0mm)
savefig(joinpath(fig_dir,"data_spectrarange_rel.png"))



## 2 - apply multiwindow spectra ##################################################

# this does a transpose as well as reduce size
data_input = transpose(data[i0_spectra:i1_spectra,:]);

# plot data subset for sanity check, compare to pump - can reuse tm_abs_spectra from above
ts_rel_spectra = map(x->( (x.instant.value)/10^9 - ts_abs_spectra[1].instant.value/10^9), ts_abs_spectra)
p_dataspectra = heatmap(ts_rel_spectra,md_values .* 3.28,data_input[:,1:di:end],fill=(true,cmap_data),clims=(-clip,clip),ylabel="Measured depth (feet)",xlabel="Time ($(time_units))",yflip=true,right_margin=5mm,colorbar=:none,xlimits=(0,t1_spectra-t0_spectra));
p_pump = plot(pumpcurve_df_forstage[:,:reldatetime],pumpcurve_df_forstage[:,"Treating Pressure"],label="",ylabel="P",right_margin=5mm,xlimits=(0,t1_spectra-t0_spectra));
l = @layout [a{0.2h} ; b]
p = plot(p_pump, p_dataspectra, layout=l)
savefig(joinpath(fig_dir,"data_spectrarange_pressure.png"))


# now do spectra

data_spectra = rfft(data_input,2) # all channels in one pass


f_nyquist = 1/(2*dt)
d_freq = f_nyquist/(size(data_spectra,2)-1)
f_vals = [(i-1)*d_freq for i in 1:size(data_spectra,2)]
f_vals_staggered = vcat([(i-1)*d_freq - d_freq/2 for i in 1:size(data_spectra,2)], [f_nyquist+d_freq/2])
data_spectra_abs = abs.(data_spectra)
di_freq = 20
spec_std = std(data_spectra_abs)
if do_bandpass
    flim= f_cut_upper
else
    flim=f_cut
end

p_spectra = heatmap(f_vals[1:di_freq:end], md_values./3.28084 .- md_values[1]/3.28084, (data_spectra_abs[:,1:di_freq:end]), xlabel="Frequency (Hz)", ylabel="Distance ($(md_units))",climits=(0,2*spec_std),fill=(true,cmap_spectra),xlimits=(0,f_cut),yflip=true)
savefig(joinpath(fig_dir,"spectra_all.png"))

p_spectra = heatmap(f_vals[1:di_freq:end], md_values./3.28084 .- md_values[1]/3.28084, log.(data_spectra_abs[:,1:di_freq:end]), xlabel="Frequency (Hz)", ylabel="Distance ($(md_units))",fill=(true,cmap_spectra),xlimits=(0,f_cut),yflip=true)
savefig(joinpath(fig_dir,"spectra_log_all.png"))

clip = std(data_spectra_abs[:,1:di_freq:end].^2)*2
p_spectra = heatmap(f_vals[1:di_freq:end], md_values./3.28084 .- md_values[1]/3.28084, (data_spectra_abs[:,1:di_freq:end].^2), xlabel="Frequency (Hz)", ylabel="Distance ($(md_units))",fill=(true,cmap_spectra),xlimits=(0,f_cut),yflip=true,clims=(0,clip))
savefig(joinpath(fig_dir,"spectra_power_all.png"))

p_spectra = heatmap(f_vals[1:di_freq:end], md_values./3.28084 .- md_values[1]/3.28084, log.(data_spectra_abs[:,1:di_freq:end].^2), xlabel="Frequency (Hz)", ylabel="Distance ($(md_units))",fill=(true,cmap_spectra),xlimits=(0,f_cut),yflip=true,clims=(0,10))
savefig(joinpath(fig_dir,"spectra_power_log_all.png"))


## 3 - extract max amplitude mode display ########################################

@info "generating spectral signal estimates for time windows..."

n_per_window = round(Int,dτ_window/dt)
n_windows = floor(Int,size(data_input,2)/n_per_window)

spec_measure = Array{Float32,2}(undef,size(data_input,1),n_windows)

for i in 1:n_windows
    # i = 1
    @info "i=$(i) of $(n_windows)"
    idx0 = (i-1)*n_per_window + 1
    idx1 = idx0 + n_per_window 
    data_freq = rfft(data_input[:,idx0:idx1],2);
    data_freq_abs = abs.(data_freq)

    n_freq = size(data_freq,2)
    f_nyquist = 1/(2dt)
    d_freq = f_nyquist/(n_freq-1)

    idx_freq0 = round(Int,f_output_min/d_freq)
    idx_freq1 = round(Int,f_output_max/d_freq)

    # signal = [ (sum(data_freq_abs[i,idx_freq0:idx_freq1] .* data_freq_abs[i,idx_freq0:idx_freq1])) for i in 1:size(data_freq_abs,1)]
    signal = [ (sum(data_freq_abs[i,idx_freq0:idx_freq1] .* data_freq_abs[i,idx_freq0:idx_freq1])^0.5) for i in 1:size(data_freq_abs,1)]
    signal = SignalUtilities.deglitch_data(signal,median_filter_length,threshold_factor=median_threshold)

    spec_measure[:,i] = signal
end

# relative times for data : time_vals[i0_spectra:di:i1_spectra]
t_spectra = [time_vals[i0_spectra] + (i-1)*dτ_window for i in 1:n_windows+1]
md_shift_stag = md_values_staggered ./ 3.28084
md_shift_stag .-= md_shift_stag[1]
md_shift_unstag = [(md_shift_stag[i] + md_shift_stag[i+1])/2 for i in 1:length(md_shift_stag)-1]
spectra_std = std(spec_measure)
minspec,maxspec=extrema(spec_measure)

heatmap(t_spectra,md_shift_stag,spec_measure,yflip=true,clims=(minspec,minspec+5*spectra_std),
        xlimits=(t_spectra[1],t_spectra[end]),fill=(true,:grays),
        xlabel="Time $(time_units)",ylabel="Distance $(md_units)")
savefig(joinpath(fig_dir,"power_$(round(Int,f_output_min))_$(round(Int,f_output_max)).png"))



stage_labels, measured_signals, md_values_dmd, tm_values_dmd, tm_abs_values_dmd = DMD.read_measures_mds_times(input_dir, idx_level);


idx_stage = 1

# find index in data read from disk by read_measures_mds_times, as that data will have all stages computed
idx_stage_raw = findfirst(x->x==get_stages[idx_stage],stage_labels)

# if findfirst(x->x==stage_labels[idx_stage], get_stages) != nothing
idx_pre_proppant = findfirst(x->x>=pre_proppant_times[idx_stage], tm_values_dmd[idx_stage_raw])

if (abs(tm_values_dmd[idx_stage_raw][idx_pre_proppant-1]-pre_proppant_times[idx_stage])
    < abs(tm_values_dmd[idx_stage_raw][idx_pre_proppant]-pre_proppant_times[idx_stage]) )
    idx_pre_proppant -= 1
end

mds_unstaggered_dmd = [(md_values_dmd[idx_stage_raw][i]+md_values_dmd[idx_stage_raw][i+1])/2 for i in 1:length(md_values_dmd[idx_stage_raw])-1]
mds_unstaggered_dmd .*= 1/3.28084
if zero_ref_mds
    md_shift = mds_unstaggered_dmd[1]
    mds_unstaggered_dmd .-= md_shift
end
if zero_ref_mds
    dist_label = "Distance ($(md_units))"
else
    dist_label = "Measured depth ($(md_units))"
end
p = plot(mds_unstaggered_dmd, measured_signals[idx_stage_raw][:,idx_pre_proppant],label="AEI",linestyle=:solid,
         xlabel=dist_label,ylabel="Acoustic signal",color=:black, xlimits=(mode_xlimits[1],mode_xlimits[2]),linewidth=2)

# now do spectra version
idx_pre_proppant = findfirst(x->x>=pre_proppant_times[idx_stage], t_spectra)
if (abs(t_spectra[idx_pre_proppant-1]-pre_proppant_times[idx_stage])
    < abs(t_spectra[idx_pre_proppant]-pre_proppant_times[idx_stage]) )
    idx_pre_proppant -= 1
end
spec_measure_normalized = spec_measure[:,idx_pre_proppant] ./ maximum(spec_measure[:,idx_pre_proppant])
plot!(p, md_shift_unstag, spec_measure_normalized, label="spectral",linestyle=:dot,color=:black)

savefig(p,joinpath(fig_dir,"compare_results_preproppant.png"))



# find index in data read from disk by read_measures_mds_times, as that data will have all stages computed

# if findfirst(x->x==stage_labels[idx_stage], get_stages) != nothing
idx_post_proppant = findfirst(x->x>=post_proppant_times[idx_stage], tm_values_dmd[idx_stage_raw])

if (abs(tm_values_dmd[idx_stage_raw][idx_post_proppant-1]-post_proppant_times[idx_stage])
    < abs(tm_values_dmd[idx_stage_raw][idx_post_proppant]-post_proppant_times[idx_stage]) )
    idx_post_proppant -= 1
end

mds_unstaggered_dmd = [(md_values_dmd[idx_stage_raw][i]+md_values_dmd[idx_stage_raw][i+1])/2 for i in 1:length(md_values_dmd[idx_stage_raw])-1]
mds_unstaggered_dmd .*= 1/3.28084
if zero_ref_mds
    md_shift = mds_unstaggered_dmd[1]
    mds_unstaggered_dmd .-= md_shift
end
if zero_ref_mds
    dist_label = "Distance ($(md_units))"
else
    dist_label = "Measured depth ($(md_units))"
end
p = plot(mds_unstaggered_dmd, measured_signals[idx_stage_raw][:,idx_post_proppant],label="AEI",linestyle=:solid,
         xlabel=dist_label,ylabel="Acoustic signal",color=:black, xlimits=(mode_xlimits[1],mode_xlimits[2]),linewidth=2)

# now do spectra version
idx_post_proppant = findfirst(x->x>=post_proppant_times[idx_stage], t_spectra)
if (abs(t_spectra[idx_post_proppant-1]-post_proppant_times[idx_stage])
    < abs(t_spectra[idx_post_proppant]-post_proppant_times[idx_stage]) )
    idx_post_proppant -= 1
end
spec_measure_normalized = spec_measure[:,idx_post_proppant] ./ maximum(spec_measure[:,idx_post_proppant])
plot!(p, md_shift_unstag, spec_measure_normalized, label="spectral",linestyle=:dot,color=:black)

savefig(p,joinpath(fig_dir,"compare_results_postproppant.png"))




# find out why times t0_spectra don't match times for DMD windows - wrong files?

