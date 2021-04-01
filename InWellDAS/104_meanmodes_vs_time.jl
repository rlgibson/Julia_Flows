#=

Apply mode extraction (mean/median) to time windows for a stage based on pump curve info.

=#
 

# assumes the following parameters are defined:

# n_levels = 10
# dτ = 10  # time window length over which mean modes will be generated and output

# # directory and file name parameters
# fig_dir = "005_crosscorr_meanmode_subtract_abs_norm_Fig"


# # data file name
# infiledir_list = ["../4001_strain_uniformity_BP_stg_501001/DMD_results",
#                   "../4002_strain_uniformity_BP_stg_501101/DMD_results",
#                   "../4003_strain_uniformity_BP_stg_501201/DMD_results"]
# #                  "../4007_strain_uniformity_stg_502301all/001_mergeStages10Hz_HDF5_99HN_5M_FRAC_99HN_STG_23_11_24_0133_DMD_results"]
# infileroot_list = ["mrDMD_100_200_HDF5_99HN_5M_FRAC_99HN_STG_10_11_20_2020_0006.h5_maxmode_",
#                    "mrDMD_100_200_HDF5_99HN_5M_FRAC_99HN_STG_11_11_20_2020_0559.h5_maxmode_",
#                    "mrDMD_100_200_HDF5_99HN_5M_FRAC_99HN_STG_12_11_20_2020_1316.h5_maxmode_"]
# #                   "mrDMD_100_200_HDF5_99HN_5M_FRAC_99HN_STG_23all_11_24_2020_0133.h5__maxmode_"]
# stage_labels = ["10","11","12"] #,"23"]

# mrdmd_dir = "./DMD_results"
# out_dir = "DMD_crosscorr_meanmode_subtract"
# outfileroot = "mrdmd_crosscorr_meanmode_subtract_"



@info "loading and compiling packages..."


using Revise, Plots, Statistics, Plots.PlotMeasures, JLD, PrettyTables, Dates, CSV, DataFrames
using DSP


# custom pkgs:
julia_path = joinpath(ENV["HOME"],"dev/nanoseis/julia_src")
# push!(LOAD_PATH,julia_path)
#push!(LOAD_PATH,joinpath(julia_path,"HDF5utilities"))
# push!(LOAD_PATH,joinpath(julia_path,"HDF5utilities/src"))
# push!(LOAD_PATH,joinpath(julia_path,"PumpCurveTools/src"))

using PumpCurveTools

using HDF5utilities

using DMD

#include("../9999_primary_flows/201_pumpcurve_tools.jl")

# include(joinpath(julia_path,"DMD/src/DMD.jl"))


# use first line to turn on debugging printouts, will get messy and slow though...
#ENV["JULIA_DEBUG"]=DMD
ENV["JULIA_DEBUG"]=""

#printstyled("debug environmental variable: $(ENV["JULIA_DEBUG"])\n",color=:magenta)

##  2 


if(!isdir(fig_dir))
    mkdir(fig_dir)
end
if(!isdir(out_dir))
    mkdir(out_dir)
end


function date_from_name(filename::String; company="halliburton")
    if company == "halliburton"
        year_range = findfirst("_2020",filename)
        year = parse(Int,filename[year_range[1]+1:year_range[end]])
        day = parse(Int,filename[year_range[1]-2:year_range[1]-1])
        month = parse(Int,filename[year_range[1]-5:year_range[1]-4])
        d = Date(year,month,day)
        @debug "date_from_name: got date $(d) from +$(filename)+"
    elseif company == "schlumberger"
        time_str = split(splitdir(filename)[2],"_")[8]
        Date(parse(Int,time_str[1:4]), parse(Int,time_str[5:6]), parse(Int,time_str[7:8]))
    else
        error("data_from_name: no file name date parser defined for company = $(company)")
    end
    return d
end

if !@isdefined company
    company = "halliburton"
end
if !@isdefined expt_date_list
    expt_date_list = map(x->date_from_name(x,company=company), infileroot_list)
end

#printstyled("date list:\n$(expt_date_list)\n",color=:brown)

if !@isdefined pumpcurve_file
    pumpcurve_df = DataFrame()
elseif pumpcurve_file == "none"
    pumpcurve_df = DataFrame()
else
    pumpcurve_df = load_pumpcurve_data(pumpcurve_file,time_shift=pumpcurve_time_shift, company=company)
    @info "***** 1st pump data time=$(pumpcurve_df[1,:truedatetime])"
end

@info "Pump curve data description:"
println(describe(pumpcurve_df))


if !@isdefined treat_pressure_hdr
    treat_pressure_hdr = "Treating Pressure"
end
if !@isdefined slurry_hdr
    slurry_hdr = "Slurry Rate"
end
if !@isdefined bh_prop_hdr
    bh_prop_hdr = "BH Proppant Conc"
end
if !@isdefined slurry_prop_hdr
    slurry_prop_hdr = "Slurry Proppant Conc"
end


@info "Starting extraction of modes..."
@info "\tbh proppant curve header=++$(bh_prop_hdr)++"
DMD.extract_modesignals_timewindows(dτ_modeextraction, ndmd_levels, out_dir, fig_dir, infileroot_list,
                                    infiledir_list, stage_labels, measure_choice, verbose_diagnostics=verbose,
                                    clip_uniformity_measure=clip_uniformity_measure,reference_scale=reference_scale,
                                    clim_scalar=clim_scalar,
                                    do_log=do_log,
                                    expt_date_list=expt_date_list, pumpcurve_df=pumpcurve_df,
                                    end_channels_omit=end_channels_omit,cmap_heatmap=cmap_dmd,
                                    md_units=md_units, md_scalar=md_scalar, zero_ref_mds=zero_ref_mds,
                                    treat_pressure_hdr = treat_pressure_hdr,
                                    slurry_hdr = slurry_hdr,
                                    bh_prop_hdr = bh_prop_hdr,
                                    slurry_prop_hdr = slurry_prop_hdr)
