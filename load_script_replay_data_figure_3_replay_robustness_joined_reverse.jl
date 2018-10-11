using HDF5

using IJulia
using PyCall
using PyPlot

@pyimport matplotlib.patches as mpatches

using JLD

d = load("data/replay_data_figure_3_with_dendritic_spikes_1.jld");


times_fast = d["times"];
ns_fast = d["ns"];
T_fast = d["T"];
ex_spike_count_per_step_per_pattern_fast = d["ex_spike_count_per_step_per_pattern"];
in_spike_count_per_step_per_pattern_fast = d["in_spike_count_per_step_per_pattern"];

d = load("data/reversal_replay_data_figure_3_with_ach_1.jld");
times_reverse = d["times"];
ns_reverse = d["ns"];
T_reverse = d["T"];
ex_spike_count_per_step_per_pattern_reverse = d["ex_spike_count_per_step_per_pattern"];
in_spike_count_per_step_per_pattern_reverse = d["in_spike_count_per_step_per_pattern"];



number_patterns = 6
Ne = 200
Ni = 100

#####smoothing out the spike count -replay
smoothing_time = 1 #divisor of T
ex_rate_over_time_smoothed_fast = zeros(number_patterns, T_fast)
in_rate_over_time_smoothed_fast = zeros(number_patterns+1, T_fast)
for time_index = 0:(T_fast-1)
	if (time_index <= smoothing_time)
		aux_spike_count = sum(ex_spike_count_per_step_per_pattern_fast[:, 1:(10*time_index+10*smoothing_time)],2) # the 10 is 1/dt of the simulation
		ex_rate_over_time_smoothed_fast[:,time_index + 1] = 1000*aux_spike_count/(Ne*(time_index+smoothing_time+1))
		aux_spike_count = sum(in_spike_count_per_step_per_pattern_fast[:, 1:(10*time_index+10*smoothing_time)],2) # the 10 is 1/dt of the simulation
		in_rate_over_time_smoothed_fast[:,time_index + 1] = 1000*aux_spike_count/(Ni*(time_index+smoothing_time+1))
	elseif time_index < T_fast - smoothing_time
		aux_spike_count = sum(ex_spike_count_per_step_per_pattern_fast[:, (10*time_index-10*smoothing_time):(10*time_index+10*smoothing_time)],2)
		ex_rate_over_time_smoothed_fast[:,time_index + 1] = 1000*aux_spike_count/(Ne*(2*smoothing_time+1))
		aux_spike_count = sum(in_spike_count_per_step_per_pattern_fast[:, (10*time_index-10*smoothing_time):(10*time_index+10*smoothing_time)],2)
		in_rate_over_time_smoothed_fast[:,time_index + 1] = 1000*aux_spike_count/(Ni*(2*smoothing_time+1))
	else
		aux_spike_count = sum(ex_spike_count_per_step_per_pattern_fast[:, (10*time_index-10*smoothing_time):(10*T_fast)],2) # the 10 is 1/dt of the simulation
		ex_rate_over_time_smoothed_fast[:,time_index + 1] = 1000*aux_spike_count/(Ne*(T_fast-time_index+smoothing_time+1))
		aux_spike_count = sum(in_spike_count_per_step_per_pattern_fast[:, (10*time_index-10*smoothing_time):(10*T_fast)],2) # the 10 is 1/dt of the simulation
		in_rate_over_time_smoothed_fast[:,time_index + 1] = 1000*aux_spike_count/(Ni*(T_fast-time_index+smoothing_time+1))
	end
end
###################################



#####smoothing out the spike count - reverse replay
smoothing_time = 1 #divisor of T
ex_rate_over_time_smoothed_reverse = zeros(number_patterns, T_reverse)
in_rate_over_time_smoothed_reverse = zeros(number_patterns+1, T_reverse)
for time_index = 0:(T_reverse-1)
	if (time_index <= smoothing_time)
		aux_spike_count = sum(ex_spike_count_per_step_per_pattern_reverse[:, 1:(10*time_index+10*smoothing_time)],2) # the 10 is 1/dt of the simulation
		ex_rate_over_time_smoothed_reverse[:,time_index + 1] = 1000*aux_spike_count/(Ne*(time_index+smoothing_time+1))
		aux_spike_count = sum(in_spike_count_per_step_per_pattern_reverse[:, 1:(10*time_index+10*smoothing_time)],2) # the 10 is 1/dt of the simulation
		in_rate_over_time_smoothed_reverse[:,time_index + 1] = 1000*aux_spike_count/(Ni*(time_index+smoothing_time+1))
	elseif time_index < T_reverse - smoothing_time
		aux_spike_count = sum(ex_spike_count_per_step_per_pattern_reverse[:, (10*time_index-10*smoothing_time):(10*time_index+10*smoothing_time)],2)
		ex_rate_over_time_smoothed_reverse[:,time_index + 1] = 1000*aux_spike_count/(Ne*(2*smoothing_time+1))
		aux_spike_count = sum(in_spike_count_per_step_per_pattern_reverse[:, (10*time_index-10*smoothing_time):(10*time_index+10*smoothing_time)],2)
		in_rate_over_time_smoothed_reverse[:,time_index + 1] = 1000*aux_spike_count/(Ni*(2*smoothing_time+1))
	else
		aux_spike_count = sum(ex_spike_count_per_step_per_pattern_reverse[:, (10*time_index-10*smoothing_time):(10*T_reverse)],2) # the 10 is 1/dt of the simulation
		ex_rate_over_time_smoothed_reverse[:,time_index + 1] = 1000*aux_spike_count/(Ne*(T_reverse-time_index+smoothing_time+1))
		aux_spike_count = sum(in_spike_count_per_step_per_pattern_reverse[:, (10*time_index-10*smoothing_time):(10*T_reverse)],2) # the 10 is 1/dt of the simulation
		in_rate_over_time_smoothed_reverse[:,time_index + 1] = 1000*aux_spike_count/(Ni*(T_reverse-time_index+smoothing_time+1))
	end
end
###################################




datestring = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
fig = plt[:figure](figsize=(14,8))
#title("Test in turning off external inhibition with a two pattern detector to prevent both patterns being active.")
ax = gca()
fig[:add_subplot](1,1,1)  # The big subplot
# Turn off axis lines and ticks of the big subplot
ax[:spines]["top"][:set_color]("none")
ax[:spines]["bottom"][:set_color]("none")
ax[:spines]["left"][:set_color]("none")
ax[:spines]["right"][:set_color]("none")
ax[:tick_params](labelcolor="w", top="off", bottom="off", left="off", right="off")
# Set common labels
ax[:set_xlabel]("Time(ms)", fontsize=35)
ax[:set_ylabel]("Ensemble", fontsize=35)



for pattern_index = 1:number_patterns
	##excitatory raster plots
	rowcount = 1

	fig[:add_subplot](number_patterns + 1, 1, pattern_index, sharex = ax) #for paper figure 3 with/without dendritic spikes
	ax1 = gca()
	setp(ax1[:get_xticklabels](),visible=false) # Disable x tick labels
	setp(ax1[:get_yticklabels](),visible=false) # Disable x tick labels
	for cc = (pattern_index-1)*Ne+1:pattern_index*Ne
		vals_fast = times_fast[cc,1:ns_fast[cc]]
		vals_reverse = times_reverse[cc,1:ns_reverse[cc]]
		y_fast = rowcount*ones(length(vals_fast))
		y_reverse = rowcount*ones(length(vals_reverse))
		scatter(vals_fast,y_fast,s=2,c="r",marker="o",linewidths=0, label="Replay")
		scatter(vals_reverse,y_reverse,s=2,c="black",marker="o",linewidths=0, label="Reverse replay")
		if pattern_index == 1 && cc == 1
			replay_patch = mpatches.Patch(color="red", label="Replay")
			reverse_replay_patch = mpatches.Patch(color="black", label="Reverse replay")
			ax1[:legend](handles=[replay_patch, reverse_replay_patch], fontsize = 15)

		end

		rowcount+=1
	end
	xlim(0,T_fast)
	ylim(0,Ne)
	if pattern_index == 1
		ylabel("\$B_1\$",fontsize=25)
	elseif pattern_index == 2
		ylabel("\$B_2\$",fontsize=25)
	elseif pattern_index == 3
		ylabel("\$B_3\$",fontsize=25)
	elseif pattern_index == 4
		ylabel("\$B_4\$",fontsize=25)
	elseif pattern_index == 5
		ylabel("\$B_5\$",fontsize=25)
	elseif pattern_index == 6
		ylabel("\$B_6\$",fontsize=25)
	end
	
end
#pattern detector inhibition raster plots

fig[:add_subplot](number_patterns + 1, 1, number_patterns + 1, sharex = ax)
ax2 = gca()
ax2[:tick_params](labelsize=20)
setp(ax2[:get_yticklabels](),visible=false) # Disable x tick labels
rowcount = 1
for cc = 1+number_patterns*Ne+number_patterns*Ni:number_patterns*Ne+(number_patterns+1)*Ni
	vals_fast = times_fast[cc,1:ns_fast[cc]]
	vals_reverse = times_reverse[cc,1:ns_reverse[cc]]	
	y_fast = rowcount*ones(length(vals_fast))
	y_reverse = rowcount*ones(length(vals_reverse))
	scatter(vals_fast,y_fast,s=2,c="b",marker="o",linewidths=0)
	scatter(vals_reverse,y_reverse,s=2,c="yellow",marker="o",linewidths=0)
	if cc == 1+number_patterns*Ne+number_patterns*Ni
		replay_patch_in = mpatches.Patch(color="blue", label="Replay")
		reverse_replay_patch_in = mpatches.Patch(color="yellow", label="Reverse replay")
		ax2[:legend](handles=[replay_patch_in, reverse_replay_patch_in], fontsize = 15)
	end
	rowcount+=1
end
ylabel("\$I_C\$",fontsize=25)
xlim(0,T_fast)
ylim(0,Ni)



tight_layout()
subplots_adjust(hspace=0.0)
savefig(string("figures/",datestring,"raster_replay_with_dendritic_spikes.png"),dpi=300)




