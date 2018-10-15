# This file is made by Marcelo Gauy and is based on
# simulation code from Litwin-Kumar 2014.

# This simulation checks whether we can turn on and off the external inhibition and get the two pattern
# detector to shut down activity in E[i] 

# Package for plotting
using PyPlot
# This is just a package to save/load from file
using HDF5

using JLD

include("simulation6.jl")
include("parameters.jl")



Ne = 200 # Number of excitatory neurons
Ni = 100 # Number of inhibitory neurons
jee0 = 10.0 # initial ee strength
jei0 = 10.0 # inhibition to excitation strength
jie = 4.0  # excitation to inhibition strength
jii = 3.0  # ii strength 
jba = 1.6 # A to B strength
jca = 1.5
jac = 12.0 
p = 0.2
rex = 1.0 # external input rate to e (khz) originally product 3
rix = 1.0 # external input rate to i (khz)
jex = 2.0  # external to e strength
jix = -4.0 # external to i strength
j_back = 0.5
ach_factor = 2.5

number_patterns = 6


function weightpars_custom() #parameters needed to generate weight matrix
	return Ne,Ni,jee0,jei0,jie,jii,jba,jca,jac,p,rex,rix,jex,jix, j_back, ach_factor
end

function stimpars_custom()
    return [1.0 100.0 0.3 1.0] #useless but the stimulus weight will only be >0 until t =100 so the synchronous stimulus must come before
end

_,_,stim_rate = stimpars_custom()

T = 400 #50 with dendritic spikes -figure 3/ 400 without figure 3
dt,_,Nskip,vpeak,Nspikes = simpars()

function simpars_custom()
    return dt, T, Nskip, vpeak, Nspikes
end

#turn_off_times = [400 800  1400  1700]
external_inhibition_turn_off_stim = [0 T 1.0 0.0]
function inhibition_1_pars()
	rextin = 1.0 # khz
	jextin = 0.0 # ???
	return rextin,jextin,external_inhibition_turn_off_stim#[turn_off_time T 1.0 0.0]
end

tic()

times,ns,times_x_axis,voltage_over_time,geaff_conductance_over_time,gerec_conductance_over_time,gi_conductance_over_time,spike_count_over_time_per_pattern,turn_on_time, ex_spike_count_per_step_per_pattern,in_spike_count_per_step_per_pattern = sim(weightpars_custom,membranepars,synapsepars,stimpars_custom,simpars_custom,inhibition_1_pars, number_patterns)
for pattern_index = 1:number_patterns
	pat_a_ex_rate = 1000*sum(ns[(pattern_index-1)*Ne+1:pattern_index*Ne])/(T*Ne)
	println("Pattern A ex rate: ",pat_a_ex_rate," Hz")
end


#####smoothing out the spike count
smoothing_time = 1 #divisor of T
ex_rate_over_time_smoothed = zeros(number_patterns, T)
in_rate_over_time_smoothed = zeros(number_patterns+1, T)
for time_index = 0:(T-1)
	if (time_index <= smoothing_time)
		aux_spike_count = sum(ex_spike_count_per_step_per_pattern[:, 1:(10*time_index+10*smoothing_time)],2) # the 10 is 1/dt of the simulation
		ex_rate_over_time_smoothed[:,time_index + 1] = 1000*aux_spike_count/(Ne*(time_index+smoothing_time+1))
		aux_spike_count = sum(in_spike_count_per_step_per_pattern[:, 1:(10*time_index+10*smoothing_time)],2) # the 10 is 1/dt of the simulation
		in_rate_over_time_smoothed[:,time_index + 1] = 1000*aux_spike_count/(Ni*(time_index+smoothing_time+1))
	elseif time_index < T - smoothing_time
		aux_spike_count = sum(ex_spike_count_per_step_per_pattern[:, (10*time_index-10*smoothing_time):(10*time_index+10*smoothing_time)],2)
		ex_rate_over_time_smoothed[:,time_index + 1] = 1000*aux_spike_count/(Ne*(2*smoothing_time+1))
		aux_spike_count = sum(in_spike_count_per_step_per_pattern[:, (10*time_index-10*smoothing_time):(10*time_index+10*smoothing_time)],2)
		in_rate_over_time_smoothed[:,time_index + 1] = 1000*aux_spike_count/(Ni*(2*smoothing_time+1))
	else
		aux_spike_count = sum(ex_spike_count_per_step_per_pattern[:, (10*time_index-10*smoothing_time):(10*T)],2) # the 10 is 1/dt of the simulation
		ex_rate_over_time_smoothed[:,time_index + 1] = 1000*aux_spike_count/(Ne*(T-time_index+smoothing_time+1))
		aux_spike_count = sum(in_spike_count_per_step_per_pattern[:, (10*time_index-10*smoothing_time):(10*T)],2) # the 10 is 1/dt of the simulation
		in_rate_over_time_smoothed[:,time_index + 1] = 1000*aux_spike_count/(Ni*(T-time_index+smoothing_time+1))
	end
end
###################################



datestring = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
fig = plt[:figure](figsize=(12,2*7))
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
	#fig[:add_subplot](2*number_patterns + 1, 1, 2*pattern_index - 1, sharex = ax)
	fig[:add_subplot](number_patterns + 1, 1, pattern_index, sharex = ax) #for paper figure 3 with/without dendritic spikes
	ax1 = gca()
	setp(ax1[:get_xticklabels](),visible=false) # Disable x tick labels
	setp(ax1[:get_yticklabels](),visible=false) # Disable x tick labels
	for cc = (pattern_index-1)*Ne+1:pattern_index*Ne
		vals = times[cc,1:ns[cc]]
		y = rowcount*ones(length(vals))
		scatter(vals,y,s=2,c="r",marker="o",linewidths=0)
		#ax[:scatter](vals,y,s=.3,c="k",marker="o",linewidths=0)
		rowcount+=1
	end
	xlim(0,T)
	ylim(0,Ne)
	if pattern_index == 1
		ylabel("\$E_1\$",fontsize=25)
	elseif pattern_index == 2
		ylabel("\$E_2\$",fontsize=25)
	elseif pattern_index == 3
		ylabel("\$E_3\$",fontsize=25)
	elseif pattern_index == 4
		ylabel("\$E_4\$",fontsize=25)
	elseif pattern_index == 5
		ylabel("\$E_5\$",fontsize=25)
	elseif pattern_index == 6
		ylabel("\$E_6\$",fontsize=25)
	end
	#ax2 = ax1[:twinx]()
	#ax2[:plot](ex_rate_over_time_smoothed[pattern_index,:],c="b")
	#ax2[:set_ylabel]("Ex. Rate", color = "b")
	#ax2[:tick_params](colors="b")
	
	#################################################
	if 0==1
	#################################################


	##inhibitory raster plots
	fig[:add_subplot](2*number_patterns + 1, 1, 2*pattern_index, sharex = ax)
	ax2 = gca()
	setp(ax2[:get_xticklabels](),visible=false) # Disable x tick labels
	setp(ax2[:get_yticklabels](),visible=false) # Disable x tick labels
	rowcount = 1
	for cc = 1+(pattern_index-1)*Ni + number_patterns*Ne:number_patterns*Ne+pattern_index*Ni
		vals = times[cc,1:ns[cc]]
		y = rowcount*ones(length(vals))
		scatter(vals,y,s=2,c="b",marker="o",linewidths=0)
		#ax[1][:scatter](vals,y,s=.3,c="k",marker="o",linewidths=0)

		rowcount+=1
	end
	xlim(0,T)
	#ylim(0,1.06*maximum(rate_over_time_smoothed[1,:]))
	ylim(0,Ni)
	if pattern_index == 1
		ylabel("\$I_1\$",fontsize=25)
	elseif pattern_index == 2
		ylabel("\$I_2\$",fontsize=25)
	elseif pattern_index == 3
		ylabel("\$I_3\$",fontsize=25)
	elseif pattern_index == 4
		ylabel("\$I_4\$",fontsize=25)
	elseif pattern_index == 5
		ylabel("\$I_5\$",fontsize=25)
	elseif pattern_index == 6
		ylabel("\$I_6\$",fontsize=25)
	end
	#ax3 = ax2[:twinx]()
	#ax3[:plot](in_rate_over_time_smoothed[pattern_index,:],c="r")
	#ax3[:set_ylabel]("In. Rate", color = "r")
	#ax3[:tick_params](colors="r")
	#xlabel("Time")
	end
end
#pattern detector inhibition raster plots

#fig[:add_subplot](2*number_patterns + 1, 1, 2*number_patterns + 1, sharex = ax)
fig[:add_subplot](number_patterns + 1, 1, number_patterns + 1, sharex = ax)
ax2 = gca()
ax2[:tick_params](labelsize=20)
setp(ax2[:get_yticklabels](),visible=false) # Disable x tick labels
rowcount = 1
for cc = 1+number_patterns*Ne+number_patterns*Ni:number_patterns*Ne+(number_patterns+1)*Ni
	vals = times[cc,1:ns[cc]]
	y = rowcount*ones(length(vals))
	scatter(vals,y,s=2,c="b",marker="o",linewidths=0)
	rowcount+=1
end
ylabel("\$I_C\$",fontsize=25)
xlim(0,T)
ylim(0,Ni)
#ax3 = ax2[:twinx]()
#ax3[:plot](in_rate_over_time_smoothed[number_patterns+1,:],c="r")
#ax3[:set_ylabel]("In. Rate", color = "r")
#ax3[:tick_params](colors="r")
#ax3[:locator_params](axis ="x", nbins=40)
#xlabel("Time")


tight_layout()
subplots_adjust(hspace=0.0)
savefig(string("figures/",datestring,"_raster_for_synchonous_with_many_pat_detector.png"),dpi=300)

toc()



save("data/replay_data_figure_3_without_dendritic_spikes.jld", "times", times, "ns", ns,"T", T, "ex_spike_count_per_step_per_pattern", ex_spike_count_per_step_per_pattern, "in_spike_count_per_step_per_pattern", in_spike_count_per_step_per_pattern)


