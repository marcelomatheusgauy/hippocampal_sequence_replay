# This file is made by Marcelo Gauy and is based on
# simulation code from Litwin-Kumar 2014.

# This simulation checks whether we can turn on and off the external inhibition and get the two pattern
# detector to shut down activity in E[i] 

# This is just a package to save/load from file
using HDF5

using IJulia
using PyCall
@pyimport matplotlib.animation as anim
using PyPlot

using JLD





include("parameters.jl")

linear_track_size = 25
loop_pos = 4
loop_length = 9
placecells_per_position = 100 
Nplacecells = placecells_per_position*linear_track_size #number of place cells
mousesize = 1
stoppingprobability = 1.0
meanwalkingdistance = 3
variancewalkingdistance = 3
minstoppingtime = 3
maxstoppingtime = 10
path_length = 100


#define the path of the mouse; 
mousepathpositions = [1.0 1.0; 2.0 1.0; 3.0 1.0; 3.0 1.0; 3.0 1.0; 3.0 1.0; 4.0 1.0; 5.0 1.0; 6.0 1.0; 7.0 1.0; 7.0 1.0; 7.0 1.0; 8.0 1.0; 9.0 1.0; 9.0 1.0; 10.0 1.0; 3.0 1.0; 4.0 1.0; 4.0 1.0; 4.0 1.0; 4.0 1.0; 4.0 1.0; 11.0 1.0; 12.0 1.0; 12.0 1.0; 12.0 1.0; 12.0 1.0; 12.0 1.0; 12.0 1.0; 12.0 1.0; 12.0 1.0; 12.0 1.0; 12.0 1.0; 12.0 1.0; 12.0 1.0; 12.0 1.0];
mousestoppingtimes = [3.0; 10.0; 14.0; 18.0; 24.0];
mousemovingtimes = [6.0; 12.0; 15.0; 22.0; 36.0]
pathlength = 36

###translate mouse path to place cell sequence

include("mousepathtocellpathsequence.jl")
cellpathsequence = placecellsequence(mousepathpositions, Nplacecells, mousesize, linear_track_size, 1)


######################################



#####################################

### get firing activity of the abstract sequence
### get place cell firing activity
### binding the sequence to the place cell activity
include("simulationsequencereplayencodingmode.jl")

Ne = 200 # Number of excitatory neurons
Ni = 100 # Number of inhibitory neurons


Ne = 200 # Number of excitatory neurons
Ni = 100 # Number of inhibitory neurons
jee0 = 10.0 # initial ee strength
jei0 = 10.0 # inhibition to excitation strength
jie = 4.0 # excitation to inhibition strength
jii = 3.0 # ii strength 
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

rex_place = 1.0
jex_place = 8.0

meanspeed = 500
variancespeed = 40

number_patterns = 25 

initial_weight_sequence_place = 0.0001

function weightpars_custom() #parameters needed to generate weight matrix
	return Ne,Ni,jee0,jei0,jie,jii,jba,jca,jac,p,rex,rix,jex,jix,j_back, Nplacecells,rex_place,jex_place,initial_weight_sequence_place
end

function stimpars_custom()
    return [1.0 24.0 0.25 30.0]
end

_,_,stim_rate = stimpars_custom()

T = 15000 ### note that this should only simulate a fraction of the path produced!



dt,_,Nskip,vpeak,Nspikes = simpars()

function simpars_custom()
    return dt, T, Nskip, vpeak, Nspikes
end

external_inhibition_turn_off_stim = T*ones(number_patterns, 4)

function inhibition_1_pars()
	rextin = 0.05 # khz
	jextin = 2.5
	return rextin,jextin,external_inhibition_turn_off_stim#[turn_off_time T 1.0 0.0]
end

#tic()

times,ns,times_x_axis,spike_count_over_time,turn_off_times,turn_on_time, ex_spike_count_per_step_per_pattern, in_spike_count_per_step_per_pattern, times_place,ns_place,time_xaxis,cellpathsequence_times, weights_to_place, times_global_inhibitory_cells, ns_global_inhibitory_cells, weights = sim(weightpars_custom,membranepars,synapsepars,stimpars_custom,simpars_custom,inhibition_1_pars, number_patterns, meanspeed,variancespeed, mousestoppingtimes,mousemovingtimes,cellpathsequence,placecells_per_position)
for pattern_index = 1:number_patterns
	pat_a_ex_rate = 1000*sum(ns[(pattern_index-1)*Ne+1:pattern_index*Ne])/(T*Ne)
	println("Pattern A ex rate: ",pat_a_ex_rate," Hz")
end
println("Encoding mode done")
println(maximum(weights_to_place))
println(minimum(weights_to_place))
#############################################
# graph for testing spiking activity of sequence
#
#####smoothing out the spike count
smoothing_time = 5 #divisor of T
ex_rate_over_time_smoothed = zeros(number_patterns, T)
in_rate_over_time_smoothed = zeros(number_patterns+1, T)
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



datestring = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
figure(figsize=(12,2*7))
title("Test in turning off external inhibition with a two pattern detector to prevent both patterns being active.")
ax = gca()
for pattern_index = 1:number_patterns
	##excitatory raster plots
	rowcount = 1
	subplot(2*number_patterns + 1, 1, 2*pattern_index - 1, sharex = ax)
	ax1 = gca()
	setp(ax1[:get_xticklabels](),visible=false) # Disable x tick labels
	for cc = (pattern_index-1)*Ne+1:pattern_index*Ne
		vals = times[cc,1:ns[cc]]
		y = rowcount*ones(length(vals))
		scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
		rowcount+=1
	end
	if number_patterns == 2
		plot([turn_off_times[1]; turn_off_times[1]],[0; T],c="r")

	else
		for index = 1:size(turn_off_times)[2]
			plot([turn_off_times[index]; turn_off_times[index]],[0; T],c="r")
			if (turn_on_time[index] != 0)
				plot([turn_on_time[index]; turn_on_time[index]],[0; T],c="b")
			end
		end
	end
	xlim(0,T)
	ylim(0,Ne)
	if pattern_index == 1
		ylabel("\$E_1\$")
	elseif pattern_index == 2
		ylabel("\$E_2\$")
	elseif pattern_index == 3
		ylabel("\$E_3\$")
	elseif pattern_index == 4
		ylabel("\$E_4\$")
	elseif pattern_index == 5
		ylabel("\$E_5\$")
	elseif pattern_index == 6
		ylabel("\$E_6\$")
	end
	ax2 = ax1[:twinx]()
	ax2[:plot](ex_rate_over_time_smoothed[pattern_index,:],c="b")
	ax2[:set_ylabel]("Ex. Rate", color = "b")
	ax2[:tick_params](colors="b")
	
	##inhibitory raster plots
	subplot(2*number_patterns + 1, 1, 2*pattern_index, sharex = ax)
	ax2 = gca()
	setp(ax2[:get_xticklabels](),visible=false) # Disable x tick labels
	rowcount = 1
	for cc = 1+(pattern_index-1)*Ni + number_patterns*Ne:number_patterns*Ne+pattern_index*Ni
		vals = times[cc,1:ns[cc]]
		y = rowcount*ones(length(vals))
		scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)

		rowcount+=1
	end
	if number_patterns == 2
		plot([turn_off_times[1]; turn_off_times[1]],[0; T],c="r")

	else
		for index = 1:size(turn_off_times)[2]
			plot([turn_off_times[index]; turn_off_times[index]],[0; T],c="r")
			if (turn_on_time[index] != 0)
				plot([turn_on_time[index]; turn_on_time[index]],[0; T],c="b")
			end
		end
	end
	xlim(0,T)

	ylim(0,Ni)
	if pattern_index == 1
		ylabel("\$I_1\$")
	elseif pattern_index == 2
		ylabel("\$I_2\$")
	elseif pattern_index == 3
		ylabel("\$I_3\$")
	elseif pattern_index == 4
		ylabel("\$I_4\$")
	elseif pattern_index == 5
		ylabel("\$I_5\$")
	elseif pattern_index == 6
		ylabel("\$I_6\$")
	end
	ax3 = ax2[:twinx]()
	ax3[:plot](in_rate_over_time_smoothed[pattern_index,:],c="r")
	in_rate_over_time_smoothed[pattern_index,:]
	ax3[:set_ylabel]("In. Rate", color = "r")
	ax3[:tick_params](colors="r")

end
#pattern detector inhibition raster plots

subplot(2*number_patterns + 1, 1, 2*number_patterns + 1, sharex = ax)
ax2 = gca()
rowcount = 1
for cc = 1+number_patterns*Ne+number_patterns*Ni:number_patterns*Ne+(number_patterns+1)*Ni
	vals = times[cc,1:ns[cc]]
	y = rowcount*ones(length(vals))
	scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
	rowcount+=1
end
if number_patterns == 2
	plot([turn_off_times[1]; turn_off_times[1]],[0; T],c="r")

else
	for index = 1:size(turn_off_times)[2]
		plot([turn_off_times[index]; turn_off_times[index]],[0; T],c="r")
		if (turn_on_time[index] != 0)
			plot([turn_on_time[index]; turn_on_time[index]],[0; T],c="b")
		end
	end
end
ylabel("\$I_C\$")
xlim(0,T)
ylim(0,Ni)
ax3 = ax2[:twinx]()
ax3[:plot](in_rate_over_time_smoothed[number_patterns+1,:],c="r")
in_rate_over_time_smoothed[number_patterns+1,:]
ax3[:set_ylabel]("In. Rate", color = "r")
ax3[:tick_params](colors="r")
xlabel("Time")


tight_layout()
subplots_adjust(hspace=0.0)
savefig(string("figures/",datestring,"_raster_for_trigger_with_many_pat_detector.png"),dpi=300)
#



#########################################################################################
### get firing activity in fast replay mode
include("simulationfastreplaymode.jl")

function weightpars_custom2() #parameters needed to generate weight matrix
	return Ne,Ni,jee0,jei0,jie,jii,jba,jca,jac,p,rex,rix,jex,jix,j_back,ach_factor, Nplacecells,1.0,0.0
end

T_fast = 150

dt,_,Nskip,vpeak,Nspikes = simpars()

function simpars_custom2()
    return dt, T_fast, Nskip, vpeak, Nspikes
end


external_inhibition_turn_off_stim = [0 T_fast 1.0 0.0]

function inhibition_1_pars2()
	rextin = 1.0 # khz
	jextin = 2.5 # ??? originally 1.4 (product)
	return rextin,jextin,external_inhibition_turn_off_stim#[turn_off_time T 1.0 0.0]
end



times_fast,ns_fast,time_xaxis_fast,spike_count_over_time_per_pattern_fast,turn_on_time_fast, ex_spike_count_per_step_per_pattern_fast,in_spike_count_per_step_per_pattern_fast,times_place_fast,ns_place_fast,time_xaxis_fast = sim(weightpars_custom2,membranepars,synapsepars,stimpars_custom,simpars_custom2,inhibition_1_pars2, number_patterns, weights_to_place, weights)

println("fast replay done")


###########################################
# produce a graphic of spikes in sequence for testing
#####smoothing out the spike count
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



datestring = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
figure(figsize=(12,2*7))
title("Test in turning off external inhibition with a two pattern detector to prevent both patterns being active.")
ax = gca()
for pattern_index = 1:number_patterns
	##excitatory raster plots
	rowcount = 1
	subplot(2*number_patterns + 1, 1, 2*pattern_index - 1, sharex = ax)
	ax1 = gca()
	setp(ax1[:get_xticklabels](),visible=false) # Disable x tick labels
	for cc = (pattern_index-1)*Ne+1:pattern_index*Ne
		vals = times_fast[cc,1:ns_fast[cc]]
		y = rowcount*ones(length(vals))
		scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)

		rowcount+=1
	end
	xlim(0,T_fast)
	ylim(0,Ne)
	if pattern_index == 1
		ylabel("\$E_1\$")
	elseif pattern_index == 2
		ylabel("\$E_2\$")
	elseif pattern_index == 3
		ylabel("\$E_3\$")
	elseif pattern_index == 4
		ylabel("\$E_4\$")
	elseif pattern_index == 5
		ylabel("\$E_5\$")
	elseif pattern_index == 6
		ylabel("\$E_6\$")
	end
	ax2 = ax1[:twinx]()
	ax2[:plot](ex_rate_over_time_smoothed_fast[pattern_index,:],c="b")
	ax2[:set_ylabel]("Ex. Rate", color = "b")
	ax2[:tick_params](colors="b")
	
	##inhibitory raster plots
	subplot(2*number_patterns + 1, 1, 2*pattern_index, sharex = ax)
	ax2 = gca()
	setp(ax2[:get_xticklabels](),visible=false) # Disable x tick labels
	rowcount = 1
	for cc = 1+(pattern_index-1)*Ni + number_patterns*Ne:number_patterns*Ne+pattern_index*Ni
		vals = times_fast[cc,1:ns_fast[cc]]
		y = rowcount*ones(length(vals))
		scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)


		rowcount+=1
	end
	xlim(0,T_fast)

	ylim(0,Ni)
	if pattern_index == 1
		ylabel("\$I_1\$")
	elseif pattern_index == 2
		ylabel("\$I_2\$")
	elseif pattern_index == 3
		ylabel("\$I_3\$")
	elseif pattern_index == 4
		ylabel("\$I_4\$")
	elseif pattern_index == 5
		ylabel("\$I_5\$")
	elseif pattern_index == 6
		ylabel("\$I_6\$")
	end
	ax3 = ax2[:twinx]()
	ax3[:plot](in_rate_over_time_smoothed_fast[pattern_index,:],c="r")
	ax3[:set_ylabel]("In. Rate", color = "r")
	ax3[:tick_params](colors="r")
	#xlabel("Time")
end
#pattern detector inhibition raster plots

subplot(2*number_patterns + 1, 1, 2*number_patterns + 1, sharex = ax)
ax2 = gca()
rowcount = 1
for cc = 1+number_patterns*Ne+number_patterns*Ni:number_patterns*Ne+(number_patterns+1)*Ni
	vals = times_fast[cc,1:ns_fast[cc]]
	y = rowcount*ones(length(vals))
	scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
	rowcount+=1
end
ylabel("\$I_C\$")
xlim(0,T_fast)
ylim(0,Ni)
ax3 = ax2[:twinx]()
ax3[:plot](in_rate_over_time_smoothed_fast[number_patterns+1,:],c="r")
ax3[:set_ylabel]("In. Rate", color = "r")
ax3[:tick_params](colors="r")
xlabel("Time")


tight_layout()
subplots_adjust(hspace=0.0)
savefig(string("figures/",datestring,"_raster_for_synchonous_with_many_pat_detector.png"),dpi=300)

#toc()


#########################################################################################
### get firing activity in fast replay mode
include("simulationreversereplay.jl")



T_reverse = 150 

dt,_,Nskip,vpeak,Nspikes = simpars()

function simpars_custom3()
    return dt, T_reverse, Nskip, vpeak, Nspikes
end



times_reverse,ns_reverse,time_xaxis_reverse,spike_count_over_time_per_pattern_reverse,turn_on_time_reverse, ex_spike_count_per_step_per_pattern_reverse,in_spike_count_per_step_per_pattern_reverse,times_place_reverse,ns_place_reverse,time_xaxis_reverse = sim(weightpars_custom2,membranepars,synapsepars,stimpars_custom,simpars_custom3,inhibition_1_pars2, number_patterns, weights_to_place, weights)

println("reverse replay done")


###########################################
# produce a graphic of spikes in sequence for testing
#####smoothing out the spike count
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
figure(figsize=(12,2*7))
title("Test in turning off external inhibition with a two pattern detector to prevent both patterns being active.")
ax = gca()
for pattern_index = 1:number_patterns
	##excitatory raster plots
	rowcount = 1
	subplot(2*number_patterns + 1, 1, 2*pattern_index - 1, sharex = ax)
	ax1 = gca()
	setp(ax1[:get_xticklabels](),visible=false) # Disable x tick labels
	for cc = (pattern_index-1)*Ne+1:pattern_index*Ne
		vals = times_reverse[cc,1:ns_reverse[cc]]
		y = rowcount*ones(length(vals))
		scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)

		rowcount+=1
	end
	xlim(0,T_reverse)
	ylim(0,Ne)
	if pattern_index == 1
		ylabel("\$E_1\$")
	elseif pattern_index == 2
		ylabel("\$E_2\$")
	elseif pattern_index == 3
		ylabel("\$E_3\$")
	elseif pattern_index == 4
		ylabel("\$E_4\$")
	elseif pattern_index == 5
		ylabel("\$E_5\$")
	elseif pattern_index == 6
		ylabel("\$E_6\$")
	end
	ax2 = ax1[:twinx]()
	ax2[:plot](ex_rate_over_time_smoothed_reverse[pattern_index,:],c="b")
	ax2[:set_ylabel]("Ex. Rate", color = "b")
	ax2[:tick_params](colors="b")
	
	##inhibitory raster plots
	subplot(2*number_patterns + 1, 1, 2*pattern_index, sharex = ax)
	ax2 = gca()
	setp(ax2[:get_xticklabels](),visible=false) # Disable x tick labels
	rowcount = 1
	for cc = 1+(pattern_index-1)*Ni + number_patterns*Ne:number_patterns*Ne+pattern_index*Ni
		vals = times_reverse[cc,1:ns_reverse[cc]]
		y = rowcount*ones(length(vals))
		scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)


		rowcount+=1
	end
	xlim(0,T_reverse)

	ylim(0,Ni)
	if pattern_index == 1
		ylabel("\$I_1\$")
	elseif pattern_index == 2
		ylabel("\$I_2\$")
	elseif pattern_index == 3
		ylabel("\$I_3\$")
	elseif pattern_index == 4
		ylabel("\$I_4\$")
	elseif pattern_index == 5
		ylabel("\$I_5\$")
	elseif pattern_index == 6
		ylabel("\$I_6\$")
	end
	ax3 = ax2[:twinx]()
	ax3[:plot](in_rate_over_time_smoothed_reverse[pattern_index,:],c="r")
	ax3[:set_ylabel]("In. Rate", color = "r")
	ax3[:tick_params](colors="r")
	#xlabel("Time")
end
#pattern detector inhibition raster plots

subplot(2*number_patterns + 1, 1, 2*number_patterns + 1, sharex = ax)
ax2 = gca()
rowcount = 1
for cc = 1+number_patterns*Ne+number_patterns*Ni:number_patterns*Ne+(number_patterns+1)*Ni
	vals = times_reverse[cc,1:ns_reverse[cc]]
	y = rowcount*ones(length(vals))
	scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
	rowcount+=1
end
ylabel("\$I_C\$")
xlim(0,T_reverse)
ylim(0,Ni)
ax3 = ax2[:twinx]()
ax3[:plot](in_rate_over_time_smoothed_reverse[number_patterns+1,:],c="r")
ax3[:set_ylabel]("In. Rate", color = "r")
ax3[:tick_params](colors="r")
xlabel("Time")


tight_layout()
subplots_adjust(hspace=0.0)
savefig(string("figures/",datestring,"_raster_for_synchonous_with_many_pat_detector.png"),dpi=300)





save("data/video_data.jld", "times", times, "ns", ns, "times_x_axis", times_x_axis, "turn_on_time", turn_on_time, "times_place", times_place,"ns_place", ns_place,"cellpathsequence_times", cellpathsequence_times,"weights_to_place",  weights_to_place, "mousepathpositions",  mousepathpositions,  "cellpathsequence",  cellpathsequence,  "times_fast",  times_fast,  "ns_fast",  ns_fast,  "times_place_fast",  times_place_fast,  "ns_place_fast",  ns_place_fast, "times_reverse",  times_reverse,  "ns_reverse",  ns_reverse,  "times_place_reverse",  times_place_reverse,  "ns_place_reverse",  ns_place_reverse, "mousestoppingtimes",  mousestoppingtimes, "mousemovingtimes",  mousemovingtimes, "linear_track_size", linear_track_size, "times_global_inhibitory_cells", times_global_inhibitory_cells, "ns_global_inhibitory_cells", ns_global_inhibitory_cells, "weights", weights)




