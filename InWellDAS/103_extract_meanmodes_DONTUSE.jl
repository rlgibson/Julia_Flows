# first pass at extracting acoustic energy indicator over time range t0_dmd to t1_dmd
# assumes the following parameters are defined:
#################################################################################
# #  time window parameters - time units are sec - recall these are relative to data in the mrDMD results, not absolutes
# t0_dmd = [2000, 2000.0, 2000] #, 3500]
# t1_dmd = [3000, 3000, 3000] #, 7500]
# n_levels = 10

# # directory and file name parameters
# fig_dir = "005_crosscorr_meanmode_subtract_abs_norm_Fig"
# clipval = 10 # max value on color scale for error heatmaps


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


# # comment this out to try to show graphics live...
# if abspath(PROGRAM_FILE) == @__FILE__
#     @info "setting turning off external display of figures (they'll only be saved to a file)"
#     ENV["GKSwstype"]=100
# end
#################################################################################


@info "loading and compiling packages..."


using Revise, Plots, Statistics, Plots.PlotMeasures, JLD, Dates, PrettyTables
using DSP

# custom pkgs:
julia_path = joinpath(ENV["HOME"],"dev/nanoseis/julia_src")
push!(LOAD_PATH,julia_path)
push!(LOAD_PATH,joinpath(julia_path,"HDF5utilities"))

using HDF5utilities


function read_mrdmd_jld(file_rootname::String,i_level::Int)
    in_levelfile = joinpath(file_rootname*"times_$(i_level).jld")
    mg_times1 = jldopen(in_levelfile,"r") do file
        read(file, "modegrid_times")
    end
    in_levelfile = joinpath(file_rootname*"times_abs_$(i_level).jld")
    mg_times_abs1 = jldopen(in_levelfile,"r") do file
        read(file, "modegrid_times_abs")
    end
    in_levelfile = joinpath(file_rootname*"chans_$(i_level).jld")
    mg_chans1 = jldopen(in_levelfile,"r") do file
        read(file, "modegrid_chans")
    end
    in_levelfile = joinpath(file_rootname*"deglitch_$(i_level).jld")
    mg_grid1 = jldopen(in_levelfile,"r") do file
        read(file, "modegrid_deglitch")
    end
    return mg_times1, mg_times_abs1, mg_chans1, mg_grid1
end


# use first line to turn on debugging printouts, will get messy and slow though...
#ENV["JULIA_DEBUG"]=Main
ENV["JULIA_DEBUG"]=""


##  2 

cd(@__DIR__)

if(!isdir(fig_dir))
    mkdir(fig_dir)
end
if(!isdir(out_dir))
    mkdir(out_dir)
end

for i_level in 1:n_levels
    @debug "i_level=$(i_level)"
     # i_level = 1
    datasets = Array{Array{Float32,2}}(undef,length(t0_dmd))
    datasets_means = Array{Array{Float32,1}}(undef,length(t0_dmd))
    datasets_ch = Array{Array{Float64,1}}(undef,length(t0_dmd))
    datasets_tm = Array{Array{Float64,1}}(undef,length(t0_dmd))
    datasets_tm_abs = Array{Array{Time,1}}(undef,length(t0_dmd))
    for i_data in 1:length(t0_dmd)
         # i_data = 1
        infileroot = joinpath(infiledir_list[i_data], infileroot_list[i_data])
        datasets_tm[i_data],datasets_tm_abs[i_data],datasets_ch[i_data],mrdmd = read_mrdmd_jld(infileroot, i_level)
        # we want to quickly find the indices of bounding times.  Since we know they're
        # ascending, cheat a little with count:
        idx_0 = length(datasets_tm[i_data])-count(>=(t0_dmd[1]), datasets_tm[i_data]) + 1
        # TODO - we may use 1 too few samples??
        idx_1 = count(<=(t1_dmd[1]), datasets_tm[i_data]) - 1
        datasets[i_data] = mrdmd[:,idx_0:idx_1]
        datasets_means[i_data] = sum(datasets[i_data],dims=2)[:,1]
        datasets_means[i_data] ./= size(datasets[i_data],2)
        # p = heatmap(tm[idx_0:idx_1+1],ch,datasets[i_data],clims=(0,3))
        # savefig(p,joinpath(fig_dir,"data_subset_level$(i_level)_group$(i_data).png"))
    end
    # TODO check selection of channel vals here
    p = plot(datasets_ch[1][2:end],datasets_means[1],label="$(stage_labels[1])",dpi=150,xlabel="Channel index",title="Mean modes, level $(i_level)")
    for i_data in 2:length(t0_dmd) #datasets_ch[i_data][2:end],
        plot!(p,datasets_ch[i_data][2:end],datasets_means[i_data],label="$(stage_labels[i_data])")
    end
    savefig(p,joinpath(fig_dir,"mean_modes_all_level$(i_level).png"))

    err_vals = Array{Float32,2}(undef,length(t0_dmd),length(t0_dmd))
    for i_data in 1:length(t0_dmd)
        for j_data in 1:length(t0_dmd)
            # i_data = j_data = 3
            @info "pair: i=$(i_data), j=$(j_data)"
            # if j_data < i_data
            #     err_vals[i_data,j_data] = err_vals[j_data,i_data]
            # else
            mag1 = sum(datasets_means[i_data] .* datasets_means[i_data])^0.5
            mag2 = sum(datasets_means[j_data] .* datasets_means[j_data])^0.5
            @debug "\tmag1=$(mag1)   mag2=$(mag2)"
            # compute cross-correlation with :longest option - this ensures 0 lag value is at
            # center of output array.  This will simplify shifting functions based on max value.
            # output length is always odd.
            xcorr_vals = DSP.xcorr(datasets_means[i_data],datasets_means[j_data]) ./ (mag1*mag2)
            @debug "\txcorr_vals extrema: $(extrema(xcorr_vals))"
            # p = plot(xcorr_vals,title="crosscorrelation: data means $(i_data) and $(j_data), level=$(i_level)")
            # savefig(p, joinpath(fig_dir,"crosscorrelation_$(i_data)_$(j_data)_level$(i_level).png"))
            
            cc_max,idx_cc_max = findmax(xcorr_vals)
            # we know xcorr result has odd length:
            idx_zerolag = round(Int,(length(xcorr_vals)-1)/2) + 1

            # if max value is to the right of zero lag, need to shift 2nd vector to the left.
            # in this case, b_shift < 0 (and vice versa)
            b_shift = idx_cc_max - idx_zerolag
            @debug "\tlength xcorr=$(length(xcorr_vals))  idx_zerolag=$(idx_zerolag)   b_shift=$(b_shift)"
            if b_shift <= 0
                n_vals = min(length(datasets_means[j_data])+b_shift, length(datasets_means[i_data]) )
                @debug "\tb_shift<=0  n_vals=$(n_vals)  b range: $(-b_shift+1):$(-b_shift+n_vals)    a range: $(1):$(n_vals)"
                err_val = (datasets_means[j_data][-b_shift+1:-b_shift+n_vals] .- datasets_means[i_data][1:n_vals]) ./ ( (mag1+mag2)/2 )
                p=plot(datasets_means[i_data][1:n_vals],label="stg=$(stage_labels[i_data])",title="Compare aligned stage mean modes")
                plot!(p,datasets_means[j_data][-b_shift+1:-b_shift+n_vals],label="stg=$(stage_labels[j_data])")
                plot!(p,err_val,label="difference")
            else
                n_vals = min(length(datasets_means[i_data])-b_shift, length(datasets_means[j_data]) )
                err_val = (datasets_means[j_data][1:n_vals] .- datasets_means[i_data][b_shift+1:b_shift+n_vals]) ./ ( (mag1+mag2)/2 )
                @debug "\tb_shift>0  n_vals=$(n_vals)    b range: $(1):$(n_vals)    a range: $(b_shift+1):$(b_shift+n_vals)"
                p=plot(datasets_means[i_data][b_shift+1:b_shift+n_vals],label="stg=$(stage_labels[i_data])",title="Compare aligned stage mean modes")
                plot!(p,datasets_means[j_data][1:n_vals],label="stg=$(stage_labels[j_data])")
                plot!(p,err_val,label="difference")
            end
            @info "\t\tnormalize by $(mag1)"
            savefig(p,joinpath(fig_dir,"meanmode_comparisons_$(i_data)_$(j_data)_level$(i_level).png"))
            # get L1 val
            err_val =  sum(abs.(err_val))
            err_vals[i_data,j_data] = err_val
            # end
        end
    end

    printstyled("Error measurements, level $(i_level):\n",color=:green)
    hl_big = Highlighter((data,i,j)->(data[i,j]>0.9),crayon"green bold")
    hl_small = Highlighter((data,i,j)->(data[i,j]<0.5),crayon"red bold")
    pretty_table(err_vals, highlighters=(hl_big,hl_small),columns_width=10)
#    formatters = ft_printf("%16.2f",1:4)
    p = heatmap(stage_labels,stage_labels,err_vals,title="uniformity measure, level $(i_level)",size=(300,300),fillcolor=:lightrainbow,xlabel="Stage",ylabel="Stage")
    savefig(p,joinpath(fig_dir,"error_matrix_level$(i_level).png"))
end
