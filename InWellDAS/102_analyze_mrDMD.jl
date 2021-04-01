#=

Read in previously computed mrDMD results.

For each level: 
1) compute mean mode and write to file, including associated md values.

2) save plots for each mean mode vs MD

3) compute uniformitymeasures for each pair of stages, using cross correlation 
   to align modes to best matching shift, subtract and normalize by average 
   mode amplitude

4) save plots of mode pairs after alignment

5) save plots of uniformity matrices for each level


------- >>  Calling code file should define the following:

>> if true, print out possibly useful stuff
verbose = false 

>> time window parameter arrays - time units are sec - recall these are relative to data in the mrDMD results, not absolutes.
>> t0 is start of time window, t1 is end
t0_dmd = [3000, # stage 1
          4000, # stage2
          ...
          2400 # stage n
          ]
t1_dmd = [5000, # stage 1
          6000, # stage 2
          ...
          4400 # stage n
          ]
n_levels = 6


>> directory and file name parameters
fig_dir = "001_crosscorr_meanmode_subtract_abs_norm_Fig"


>> input data directory list (array)
infiledir_list = ["../3001_strain_uniformity_stg_501001/DMD_results", # stage 1
                  "../3002_strain_uniformity_stg_501101/DMD_results", # stage 2
                  ...
                  "../3003o_strain_uniformity_stg_503601/DMD_results" # stage n
                  ]
>> root file names from mrDMD work 
infileroot_list = ["mrDMD_f10Hz_HDF5_99HN_5M_FRAC_99HN_STG_10_11_20_2020_0006.h5_maxmode_", # stage 1
                   "mrDMD_f10Hz_HDF5_99HN_5M_FRAC_99HN_STG_11_11_20_2020_0559.h5_maxmode_", # stage 2 
                   ...
                   "mrDMD_f10Hz_HDF5_99HN_5M_FRAC_99HN_STG_36_11_29_2020.h5_maxmode_" # stage n
                   ]
>> labels (names) for each stage
stage_labels = ["1", "2", ..., "n"]

>> computed output files saved in this directory
out_dir = "DMD_crosscorr_meanmode_subtract"

=#

## 1

@info "loading and compiling packages..."


using Revise, Plots, Statistics, Plots.PlotMeasures, JLD, PrettyTables, Dates
using DSP

# custom pkgs:
julia_path = joinpath(ENV["HOME"],"dev/nanoseis/julia_src")
push!(LOAD_PATH,julia_path)
push!(LOAD_PATH,joinpath(julia_path,"HDF5utilities/src"))

@info LOAD_PATH

include(joinpath(julia_path,"DMD/src/DMD.jl"))

using HDF5utilities

#include("../9999_primary_flows/200_mrDMD_analysis_tools.jl")

# use first line to turn on debugging printouts, will get messy and slow though...
# ENV["JULIA_DEBUG"]=Main
# ENV["JULIA_DEBUG"]="DMD"


##  2 

if length(t0_dmd) != length(t1_dmd) #
    error("Mismatch in number of data values in t0_dmd and t1_dmd")
end
if  length(t0_dmd) != length(infiledir_list)
    error("Mismatch in number of data values in t0_dmd and infiledir_list")
end
if length(t0_dmd) != length(infileroot_list)
    error("Mismatch in number of data values in t0_dmd and infileroot_list")
end
if  length(t0_dmd) != length(stage_labels)
    error("Mismatch in number of data values in t0_dmd and stage_labels")
end

# cd(@__DIR__)

if(!isdir(fig_dir))
    mkdir(fig_dir)
end
if(!isdir(out_dir))
    mkdir(out_dir)
end

@info "processing mrDMD..."

extracted_modes, datasets_ch, datasets_tm, datasets_tm_abs, datasets_md_values, mrdmd_results, md_values_unstaggered =
    DMD.extract_modesignals(t0_dmd, t1_dmd, n_levels, out_dir, fig_dir, infileroot_list,
                            infiledir_list, stage_labels, measure_choice, verbose_diagnostics=verbose,
                            clip_uniformity_measure=clip_uniformity_measure)

