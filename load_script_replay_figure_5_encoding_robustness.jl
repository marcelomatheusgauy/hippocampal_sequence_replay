using HDF5

using IJulia
using PyCall
using PyPlot

using JLD

d = load("data/replay_data_1.jld");

times = d["times"];
ns = d["ns"];
times_x_axis = d["times_x_axis"];
turn_on_time = d["turn_on_time"];
turn_off_times = d["turn_off_times"];
ex_spike_count_per_step_per_pattern = d["ex_spike_count_per_step_per_pattern"];
in_spike_count_per_step_per_pattern = d["in_spike_count_per_step_per_pattern"];
times_global_inhibitory_cells = d["times_global_inhibitory_cells"];
ns_global_inhibitory_cells = d["ns_global_inhibitory_cells"];


number_patterns = 6
Ne = 200
Ni = 100

global_inhibitory_cells = 200

Nsteps = size(times_x_axis,1)
T = convert(Int64,Nsteps*0.1)

println(Nsteps)
println(T)






#####smoothing out the spike count
smoothing_time = 10 #divisor of T
ex_rate_over_time_smoothed = zeros(number_patterns, T)
in_rate_over_time_smoothed = zeros(number_patterns+2, T)
for time_index = 1:T
	if (time_index <= smoothing_time)
		aux_spike_count = sum(ex_spike_count_per_step_per_pattern[:, 1:(10*time_index+10*smoothing_time)],2) # the 10 is 1/dt of the simulation
		ex_rate_over_time_smoothed[:,time_index] = 1000*aux_spike_count/(Ne*(time_index+smoothing_time+1))
		aux_spike_count = sum(in_spike_count_per_step_per_pattern[:, 1:(10*time_index+10*smoothing_time)],2) # the 10 is 1/dt of the simulation
		in_rate_over_time_smoothed[:,time_index] = 1000*aux_spike_count/(Ni*(time_index+smoothing_time+1))
	elseif time_index < T - smoothing_time
		aux_spike_count = sum(ex_spike_count_per_step_per_pattern[:, (10*time_index-10*smoothing_time):(10*time_index+10*smoothing_time)],2)
		ex_rate_over_time_smoothed[:,time_index] = 1000*aux_spike_count/(Ne*(2*smoothing_time+1))
		aux_spike_count = sum(in_spike_count_per_step_per_pattern[:, (10*time_index-10*smoothing_time):(10*time_index+10*smoothing_time)],2)
		in_rate_over_time_smoothed[:,time_index] = 1000*aux_spike_count/(Ni*(2*smoothing_time+1))
	else
		aux_spike_count = sum(ex_spike_count_per_step_per_pattern[:, (10*time_index-10*smoothing_time):(10*T)],2) # the 10 is 1/dt of the simulation
		ex_rate_over_time_smoothed[:,time_index] = 1000*aux_spike_count/(Ne*(T-time_index+smoothing_time+1))
		aux_spike_count = sum(in_spike_count_per_step_per_pattern[:, (10*time_index-10*smoothing_time):(10*T)],2) # the 10 is 1/dt of the simulation
		in_rate_over_time_smoothed[:,time_index] = 1000*aux_spike_count/(Ni*(T-time_index+smoothing_time+1))
	end
end
###################################

#####spike count
spike_count_global_inhibitory_cells_per_time = zeros(T)
spike_count_global_inhibitory_cells = zeros(Int64, global_inhibitory_cells)
for time_index = 1:T
	for global_inhibitory_index in 1:global_inhibitory_cells
		cc = global_inhibitory_index
		next_spike = spike_count_global_inhibitory_cells[cc] + 1
		if next_spike < ns_global_inhibitory_cells[cc]
			if 0.0 < times_global_inhibitory_cells[cc,next_spike+1] <= time_index
				spike_count_global_inhibitory_cells_per_time[time_index] += 1.0
				spike_count_global_inhibitory_cells[cc] += 1
			end
		end
	end
end



###################################
#####smoothing out the spike count
for time_index = 1:T
	if time_index <= smoothing_time
		aux_spike_count = sum(spike_count_global_inhibitory_cells_per_time[1:(time_index+smoothing_time)],1) # the 10 is 1/dt of the simulation
		in_rate_over_time_smoothed[number_patterns+2,time_index] = 1000*aux_spike_count[1]/(global_inhibitory_cells*(time_index+smoothing_time+1))
	elseif time_index < T - smoothing_time
		aux_spike_count = sum(spike_count_global_inhibitory_cells_per_time[(time_index-smoothing_time):(time_index+smoothing_time)],1)
		in_rate_over_time_smoothed[number_patterns+2,time_index] = 1000*aux_spike_count[1]/(global_inhibitory_cells*(2*smoothing_time+1))
	else
		aux_spike_count = sum(spike_count_global_inhibitory_cells_per_time[(time_index-smoothing_time):T],1) # the 10 is 1/dt of the simulation
		in_rate_over_time_smoothed[number_patterns+2,time_index]=1000*aux_spike_count[1]/(global_inhibitory_cells*(T-time_index+smoothing_time+1))
	end
end

####################################


datestring = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
fig = plt[:figure](figsize=(16,14))
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
ax[:set_xlabel]("Time(ms)", fontsize=40)
ax[:set_ylabel]("Ensemble", fontsize=40)
ax[:yaxis][:labelpad] = 55
axsec = ax[:twinx]()
axsec[:set_ylabel]("Rate(Hz)", fontsize=40)
axsec[:yaxis][:labelpad] = 30
axsec[:spines]["top"][:set_color]("none")
axsec[:spines]["bottom"][:set_color]("none")
axsec[:spines]["left"][:set_color]("none")
axsec[:spines]["right"][:set_color]("none")
axsec[:tick_params](labelcolor="w", top="off", bottom="off", left="off", right="off")

for pattern_index = 1:number_patterns
	##excitatory raster plots
	rowcount = 1

	fig[:add_subplot](number_patterns + 2, 1, pattern_index, sharex = ax) #2 patterns figure 2D without inhibitory ensembles
	ax1 = gca()
	setp(ax1[:get_xticklabels](),visible=false) # Disable x tick labels
	setp(ax1[:get_yticklabels](),fontsize=20) # Disable x tick labels

	for cc = (pattern_index-1)*Ne+1:pattern_index*Ne
		vals = times[cc,1:ns[cc]]
		y = rowcount*ones(length(vals))
		scatter(vals,y,s=2,c="k",marker="o",linewidths=0)

		rowcount+=1
	end
	xlim(0,T)
	ylim(0,Ne)
	if pattern_index == 1
		ylabel("\$B_1\$",fontsize=30) # figure 2D 2 patterns no inhibition in ensembles
	elseif pattern_index == 2
		ylabel("\$B_2\$",fontsize=30) # figure 2D 2 patterns no inhibition in ensembles
	elseif pattern_index == 3
		ylabel("\$B_3\$",fontsize=30)
	elseif pattern_index == 4
		ylabel("\$B_4\$",fontsize=30)
	elseif pattern_index == 5
		ylabel("\$B_5\$",fontsize=30)
	elseif pattern_index == 6
		ylabel("\$B_6\$",fontsize=30)
	end
	ax2 = ax1[:twinx]()
	ax2[:plot](ex_rate_over_time_smoothed[pattern_index,:],c="r")
	ax2[:tick_params](colors="r")

	ax2[:tick_params](labelsize = 20)
	ax2[:tick_params](axis="y", which="both", left="off")
	ylim(0,80)
	xlim(0,T)
end
#pattern detector inhibition raster plots


fig[:add_subplot](number_patterns + 2, 1, number_patterns + 1, sharex = ax) # 2 patterns figure 2D without inhibitory ensembles
ax2 = gca()
setp(ax2[:get_xticklabels](),visible=false) # Disable x tick labels
setp(ax2[:get_yticklabels](),fontsize=20) # Disable x tick labels

rowcount = 1
for cc = 1+number_patterns*Ne+number_patterns*Ni:number_patterns*Ne+(number_patterns+1)*Ni
	vals = times[cc,1:ns[cc]]
	y = rowcount*ones(length(vals))
	scatter(vals,y,s=2,c="c",marker="o",linewidths=0)
	rowcount+=1
end
ylabel("\$I_C\$",fontsize=30)

ylim(0,Ni)
ax3 = ax2[:twinx]()
ax3[:plot](in_rate_over_time_smoothed[number_patterns+1,:],c="c")


ax3[:tick_params](colors="c")
ax3[:tick_params](labelsize = 20)



fig[:add_subplot](number_patterns + 2, 1, number_patterns + 2, sharex = ax) # 2 patterns figure 2D without inhibitory ensembles
ax4 = gca()
ax4[:tick_params](labelsize=20)

setp(ax4[:get_yticklabels](),fontsize = 20) # Disable x tick labels
rowcount = 1
for cc = 1:200
	vals = times_global_inhibitory_cells[cc,1:ns_global_inhibitory_cells[cc]]
	y = rowcount*ones(length(vals))
	scatter(vals,y,s=2,c="b",marker="o",linewidths=0)
	rowcount+=1
end
ylabel("\$I_G\$",fontsize=30)

ylim(0,global_inhibitory_cells)
ax5 = ax4[:twinx]()
ax5[:plot](in_rate_over_time_smoothed[number_patterns + 2,:],c="b")

ax5[:tick_params](colors="b")
ax5[:tick_params](labelsize = 20)
ax5[:tick_params](axis="y", which="both", left="off")




tight_layout()
subplots_adjust(hspace=0.0)
savefig(string("figures/",datestring,"spiking_two_patterns_dynamics_zoom.png"),dpi=300) #figure 2D



