using DSP, JLD, Plots

"""
    read_mrdmd_jld(file_rootname::String,i_level::Int)

# Return:
mg_times1, mg_times_abs1, mg_chans1, mg_md_stag, mg_grid1
"""
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
    in_levelfile = joinpath(file_rootname*"md_values_stagger_$(i_level).jld")
    mg_md_stag = jldopen(in_levelfile,"r") do file
        read(file, "modegrid_md_values_stagger")
    end
    return mg_times1, mg_times_abs1, mg_chans1, mg_md_stag, mg_grid1
end



function uniformity_measure(t0_dmd::Array{Float64,1}, t1_dmd::Array{Float64,1}, mode_vals::Array{Array{Float32,1},1},
                            idx_level::Int, fig_dir::String, stage_labels::Array{String,1};
                            verbose_diagnostics::Bool=false
                            )::Array{Float32,2}

    err_vals = Array{Float32,2}(undef,length(t0_dmd),length(t0_dmd))

    for i_data in 1:length(t0_dmd)
        for j_data in 1:length(t0_dmd)
            if verbose_diagnostics
                @info "pair: i=$(i_data), j=$(j_data)"
            end

            # get magnitude of each mode for normalizing results
            mag1 = sum(mode_vals[i_data] .* mode_vals[i_data])^0.5
            mag2 = sum(mode_vals[j_data] .* mode_vals[j_data])^0.5

            # compute cross-correlation with :longest option - this ensures 0 lag value is at
            # center of output array.  This will simplify shifting functions based on max value.
            # output length is always odd.
            xcorr_vals = DSP.xcorr(mode_vals[i_data],mode_vals[j_data]) ./ (mag1*mag2)

            cc_max,idx_cc_max = findmax(xcorr_vals)
            idx_zerolag = round(Int,length(xcorr_vals)/2) 

            # if max value is to the right of zero lag, need to shift 2nd vector to the left.
            # in this case, b_shift < 0 (and vice versa)
            b_shift = idx_cc_max - idx_zerolag
            if b_shift <= 0
                n_vals = min(length(mode_vals[j_data])+b_shift, length(mode_vals[i_data]) )
                err_val = (mode_vals[j_data][-b_shift+1:-b_shift+n_vals] .- mode_vals[i_data][1:n_vals]) ./ ( (mag1+mag2)/2 )
                p=plot(mode_vals[i_data][1:n_vals],label="stg=$(stage_labels[i_data])",title="Compare aligned stage mean modes",
                       xlabel="Arbitrary channel index post alignment",ylabel="Acoustic energy")
                plot!(p,mode_vals[j_data][-b_shift+1:-b_shift+n_vals],label="stg=$(stage_labels[j_data])")
                plot!(p,err_val,label="difference")
            else
                n_vals = min(length(mode_vals[i_data])-b_shift, length(mode_vals[j_data]) )
                err_val = (mode_vals[j_data][1:n_vals] .- mode_vals[i_data][b_shift+1:b_shift+n_vals]) ./ ( (mag1+mag2)/2 )
                p=plot(mode_vals[i_data][b_shift+1:b_shift+n_vals],label="stg=$(stage_labels[i_data])",title="Compare aligned stage mean modes",
                       xlabel="Arbitrary channel index post alignment",ylabel="Acoustic energy")
                plot!(p,mode_vals[j_data][1:n_vals],label="stg=$(stage_labels[j_data])")
                plot!(p,err_val,label="difference")
            end

            if verbose_diagnostics
                @info "\t\tnormalize by $((mag1+mag2)/2)"
            end

            savefig(p,joinpath(fig_dir,"meanmode_comparisons_$(i_data)_$(j_data)_level$(idx_level).png"))

            # get L1 val
            err_val =  sum(abs.(err_val))
            err_vals[i_data,j_data] = err_val
        end
    end

    return err_vals
end


function plot_all_modes(md_values_unstaggered::Array{Array{Float32,1},1}, extracted_modes::Array{Array{Float32,1},1}, idx_level::Int,
                        fig_dir::String, stage_labels::Array{String,1})

    # TODO check selection of channel vals here
    p = plot(md_values_unstaggered[1], extracted_modes[1], label="$(stage_labels[1])",dpi=150,xlabel="Measured depth (feet)",
             title="Mean modes, level $(idx_level)")
    for i_data in 2:length(extracted_modes)
        plot!(p,md_values_unstaggered[i_data], extracted_modes[i_data],label="$(stage_labels[i_data])")
    end
    savefig(p,joinpath(fig_dir,"mean_modes_all_level$(idx_level).png"))

end



function get_data_and_measures(extracted_modes::Array{Array{Float32,1},1},
                               datasets_ch::Array{Array{Float64,1},1},
                               datasets_tm::Array{Array{Float64,1},1},
                               datasets_tm_abs::Array{Array{Time,1},1},
                               datasets_md_values::Array{Array{Float32,1},1},
                               mrdmd_results::Array{Array{Float32,2},1},
                               md_values_unstaggered::Array{Array{Float32,1},1},
                               idx_level::Int,
                               infiledir_list::Array{String,1},
                               infileroot_list::Array{String,1}
)

    for i_data in 1:length(extracted_modes)

        # i_data = 1

        infileroot = joinpath(infiledir_list[i_data], infileroot_list[i_data])
        
        @debug "i_data=$(i_data)   idx_level=$(idx_level)   about to call reader for infileroot=$(infileroot)"
        datasets_tm[i_data], datasets_tm_abs[i_data], datasets_ch[i_data], datasets_md_values[i_data], mrdmd_results[i_data] =
            read_mrdmd_jld(infileroot, idx_level)

        @debug "after read, sizes: tm $(size(datasets_tm[i_data])), tm_abs $(size(datasets_tm_abs[i_data])), ch $(size(datasets_ch[i_data])), md $(size(datasets_md_values[i_data])), dmd: $(size(mrdmd))"
        md_values_unstaggered[i_data] = [(datasets_md_values[i_data][i] + datasets_md_values[i_data][i])/2 for i in 1:length(datasets_md_values[i_data])-1]
        md_values_unstaggered[i_data] # should be in ft already :-(
        
        # we want to quickly find the indices of bounding times.  Since we know they're
        # ascending, cheat a little with count:
        idx_0 = length(datasets_tm[i_data])-count(>=(t0_dmd[i_data]), datasets_tm[i_data]) + 1
        # TODO - we may use 1 too few samples??
        idx_1 = count(<=(t1_dmd[1]), datasets_tm[i_data]) - 1

        @debug "t0_dmd[i_data] = $(t0_dmd[i_data]),  1st dataset tm=$(datasets_tm[i_data][1]), last=$(datasets_tm[i_data][end])"
        # datasets[i_data] = mrdmd_results[i_data][:,idx_0:idx_1]
        if measure_choice == :mean
            extracted_modes[i_data] = mean(mrdmd_results[i_data][:,idx_0:idx_1], dims=2)[:,1]
        elseif measure_choice == :median
            extracted_modes[i_data] = median(mrdmd_results[i_data][:,idx_0:idx_1], dims=2)[:,1]
        end

        @debug "idx_level=$(idx_level), i_data=$(i_data), extrema of mrDMD result: $(extrema(mrdmd_results[i_data]))"
        @debug "\textrema of measured signal vector: $(extrema(extracted_modes[i_data]))"


    end

end



function save_measures_mds(measured_signals::Array{Array{Float32,1},1}, md_values::Array{Array{Float32,1},1},
                           idx_level::Int, out_dir::String)

    for i_data in 1:length(measured_signals)
        jldopen(joinpath(out_dir,"mean_mode_level$(idx_level)_stg$(stage_labels[i_data]).jld"),"w") do file
            write(file, "mean_mode",measured_signals[i_data])
        end
        jldopen(joinpath(out_dir,"mds_for_means_level$(idx_level)_$(i_data).jld"),"w") do file
            write(file,"mds",md_values[i_data])
        end
    end
    
end


function extract_modesignals(t0_dmd::Array{Float64,1},
                               t1_dmd::Array{Float64,1},
                               out_dir::String,
                               fig_dir::String,
                               stage_labels::Array{String,1};
                               verbose_diagnostics::Bool=false)

    extracted_modes = Array{Array{Float32,1}}(undef,length(t0_dmd))
    datasets_ch = Array{Array{Float64,1}}(undef,length(t0_dmd))
    datasets_tm = Array{Array{Float64,1}}(undef,length(t0_dmd))
    datasets_tm_abs = Array{Array{Time,1}}(undef,length(t0_dmd))
    datasets_md_values = Array{Array{Float32,1}}(undef,length(t0_dmd))
    mrdmd_results = Array{Array{Float32,2}}(undef,length(infiledir_list))
    md_values_unstaggered = Array{Array{Float32,1}}(undef,length(t0_dmd))

    for i_level in 1:n_levels

        # i_level = 1

        get_data_and_measures(extracted_modes, datasets_ch, datasets_tm, datasets_tm_abs, datasets_md_values,
                              mrdmd_results, md_values_unstaggered, i_level, infiledir_list,
                              infileroot_list)

        save_measures_mds(extracted_modes, md_values_unstaggered, i_level, out_dir)

        plot_all_modes(md_values_unstaggered, extracted_modes, i_level, fig_dir, stage_labels)

        uniformity_vals = uniformity_measure(t0_dmd, t1_dmd, extracted_modes, i_level, fig_dir, stage_labels, verbose_diagnostics=verbose_diagnostics)

        printstyled("Error measurements, level $(i_level):\n",color=:green)
        hl_big = Highlighter((data,i,j)->(data[i,j]>0.9),crayon"green bold")
        hl_small = Highlighter((data,i,j)->(data[i,j]<0.5),crayon"red bold")
        pretty_table(uniformity_vals, highlighters=(hl_big,hl_small),columns_width=10)

        p = heatmap(stage_labels,stage_labels,uniformity_vals,title="stage uniformity matrix, level $(i_level)",size=(300,300),
                    fillcolor=:lightrainbow,xlabel="Stage",ylabel="Stage")
        savefig(p,joinpath(fig_dir,"uniformity_matrix_level$(i_level).png"))
    end

    return extracted_modes, datasets_ch, datasets_tm, datasets_tm_abs, datasets_md_values, mrdmd_results, md_values_unstaggered 
end
