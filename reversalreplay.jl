# Package for plotting
using PyPlot
# This is just a package to save/load from file
using HDF5
using JLD

include("simulation8.jl")
include("parameters.jl")

# Ne = 200 # Number of excitatory neurons
# Ni = 100 # Number of inhibitory neurons
# jee0 = 6.0 # initial ee strength
# jei0 = 7.0 # inhibition to excitation strength
# jie = 5.0 # excitation to inhibition strength
# jii = 5.0 # ii strength 
# jba = 0.5 # A to B strength
# p = 0.2
# rex = 0.6 # external input rate to e (khz) 
# rix = 2.0 # external input rate to i (khz)
# jex = 6.0  # external to e strength
# jix = -6.0 # external to i strength

noise = rand(8)
noise = ones(8)


Ne = 200 # Number of excitatory neurons
Ni = 100 # Number of inhibitory neurons
jee0 = 9.8 + 0.2*noise[1] # initial ee strength
jei0 = 9.8 + 0.2*noise[2] # inhibition to excitation strength
jie = 3.8 + 0.2*noise[3] # excitation to inhibition strength
jii = 2.8 + 0.2*noise[4] # ii strength 
jba = 1.4 + 0.2*noise[5]# A to B strength
jca = 1.3 + 0.2*noise[7]
jac = 11.8 + 0.2*noise[8]
p = 0.2
rex = 1.0 # external input rate to e (khz) originally product 3
rix = 1.0 # external input rate to i (khz)
jex = 2.0  # external to e strength
jix = -4.0 # external to i strength
j_back = 0.5
ach_factor = 1.0

number_patterns = 6


function weightpars_custom() #parameters needed to generate weight matrix
	return Ne,Ni,jee0,jei0,jie,jii,jba,jca,jac,p,rex,rix,jex,jix, j_back,ach_factor
end

function stimpars_custom()
    return [1.0 100.0 0.3 20.0] #useless but the stimulus weight will only be >0 until t =100 so the synchronous stimulus must come before
end

_,_,stim_rate = stimpars_custom()

T = 50
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

#pat_a_ex_rate = 1000*sum(ns[1:Ne])/(T*Ne)
#pat_b_ex_rate = 1000*sum(ns[Ne+1:2*Ne])/(T*Ne)
#pat_c_ex_rate = 1000*sum(ns[2*Ne+1:3*Ne])/(T*Ne)
#println("Pattern A ex rate: ",pat_a_ex_rate," Hz")
#println("Pattern B ex rate: ",pat_b_ex_rate," Hz")
#println("Pattern C ex rate: ",pat_c_ex_rate," Hz")


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
#axsec = ax[:twinx]()
#axsec[:set_ylabel]("Rate(Hz)", fontsize=35)
#axsec[:spines]["top"][:set_color]("none")
#axsec[:spines]["bottom"][:set_color]("none")
#axsec[:spines]["left"][:set_color]("none")
#axsec[:spines]["right"][:set_color]("none")
#axsec[:tick_params](labelcolor="w", top="off", bottom="off", left="off", right="off")
#figure(figsize=(12,2*7))
#title("Test in turning off external inhibition with a two pattern detector to prevent both patterns being active.")
#ax = gca()
for pattern_index = 1:number_patterns
	##excitatory raster plots
	rowcount = 1
	fig[:add_subplot](2*number_patterns + 1, 1, 2*pattern_index - 1, sharex = ax)
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
	#ax2 = ax1[:twinx]()
	#ax2[:plot](ex_rate_over_time_smoothed[pattern_index,:],c="b")
	#ax2[:set_ylabel]("Ex. Rate", color = "b")
	#ax2[:tick_params](colors="b")
	
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
#pattern detector inhibition raster plots

fig[:add_subplot](2*number_patterns + 1, 1, 2*number_patterns + 1, sharex = ax)
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
savefig(string("figures/",datestring,"_raster_for_synchonous_reverse_with_many_pat_detector.png"),dpi=300)

toc()



save("data/reversal_replay_data_figure_3_with_ach.jld", "times", times, "ns", ns,"T", T, "ex_spike_count_per_step_per_pattern", ex_spike_count_per_step_per_pattern, "in_spike_count_per_step_per_pattern", in_spike_count_per_step_per_pattern)


#=
datestring = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
figure(figsize=(12,2*7))
title("Test in turning off external inhibition with a two pattern detector to prevent both patterns being active.")
rowcount = 1
for pattern_index = 1:number_patterns
	for cc = (pattern_index-1)*Ne+1:pattern_index*Ne
		subplot(2*number_patterns+1,1,pattern_index)
		vals = times[cc,1:ns[cc]]
		y = rowcount*ones(length(vals))
		scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
		rowcount+=1
	end
	#for index = 1:size(turn_off_times)[2]
	#	plot([turn_off_times[index]; turn_off_times[index]],[0; T],c="r")
	#end
	ylabel("PE")
	xlim(0,T)
	ylim(0,Ne)
	rowcount = 1
end
for pattern_index = 1:number_patterns
	for cc = 1+(pattern_index-1)*Ni + number_patterns*Ne:number_patterns*Ne+pattern_index*Ni
		subplot(2*number_patterns+1,1,number_patterns + pattern_index)
		vals = times[cc,1:ns[cc]]
		y = rowcount*ones(length(vals))
		scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
		rowcount+=1
	end
	#for index = 1:size(turn_off_times)[2]
	#	plot([turn_off_times[index]; turn_off_times[index]],[0; T],c="r")
	#end
	ylabel("PI")
	xlim(0,T)
	ylim(0,Ni)
	rowcount = 1
end

for cc = 1+number_patterns*Ne+number_patterns*Ni:number_patterns*Ne+(number_patterns+1)*Ni
	subplot(2*number_patterns+1,1,2*number_patterns+1)
	vals = times[cc,1:ns[cc]]
	y = rowcount*ones(length(vals))
	scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
	rowcount+=1
end
#for index = 1:size(turn_off_times)[2]
#	plot([turn_off_times[index]; turn_off_times[index]],[0; T],c="r")
#end
ylabel("PDI")
xlim(0,T)
ylim(0,Ni)
#xlim(0,T)
#ylim(0,2*Ne + 3*Ni)
#ylabel("Neuron")
#xlabel("Time")
tight_layout()
savefig(string("figures/",datestring,"_raster_for_synchonous_with_many_pat_detector.png"),dpi=300)

toc()



#clf()
#figure(figsize=(12,2*4))
#subplot(4,1,1)
#plot(transpose(voltage_over_time[1,:]),c="b")
#plot(transpose(voltage_over_time[3,:]),c="r")
#ylabel("Membrane voltage")
#title("Pattern A, neuron 1")
#ylim(-66,-10)

#subplot(4,1,2)
#plot(transpose(voltage_over_time[2,:]),c="b")
#ylabel("Membrane voltage")
#title("Pattern A, neuron 2")
#ylim(-66,-10)

#subplot(4,1,3)
#plot(transpose(voltage_over_time[3,:]),c="b")
#ylabel("Membrane voltage")
#title("Pattern B, neuron 1")
#ylim(-66,-10)

#subplot(4,1,4)
#plot(transpose(voltage_over_time[4,:]),c="b")
#ylabel("Membrane voltage")
#title("Pattern B, neuron 2")
#ylim(-66,-10)

#tight_layout()
#savefig(string("figures/",datestring,"_voltage_for_synchoronous.png"),dpi=300)

#println(ge_conductance_over_time[1,30])
clf()
figure(figsize=(12,2*4))
subplot(4,1,1)
plot(transpose(geaff_conductance_over_time[1,:]),c="b")
plot(transpose(geaff_conductance_over_time[Ne+1,:]),c="r")
ylabel("Excitatory aff conductance")
title("Pattern A, neuron 1")
#ylim(0,5)

subplot(4,1,2)
plot(transpose(geaff_conductance_over_time[Ne+1,:]),c="b")
ylabel("Excitatory aff conductance")
title("Pattern B, neuron 1")
#ylim(0,5)

subplot(4,1,3)
ax=gca()
plot(transpose(mean(geaff_conductance_over_time[5*Ne+1:6*Ne,:],1)),c="b")
plot(transpose(mean(geaff_conductance_over_time[4*Ne+1:5*Ne,:],1)),c="r")
ax[:locator_params](axis ="x", nbins=40)
ylabel("Excitatory aff conductance")
title("Pattern A, average")
#ylim(0,5)

subplot(4,1,4)
plot(transpose(mean(geaff_conductance_over_time[Ne+1:2*Ne,:],1)),c="b")
plot(transpose(mean(geaff_conductance_over_time[2*Ne+1:3*Ne,:],1)),c="r")
plot(transpose(mean(geaff_conductance_over_time[3*Ne+1:4*Ne,:],1)),c="g")
plot(transpose(mean(geaff_conductance_over_time[4*Ne+1:5*Ne,:],1)),c="k")
ylabel("Excitatory aff conductance")
title("Pattern B, average")

#ylim(0,5)
#subplot(6,1,5)
#plot(transpose(gi_conductance_over_time[Ne*number_patterns + 1,:]),c="b")
#ylabel("Inhibitory conductance")
#title("Pattern A, inh neuron 1")
#ylim(-66,-10)

#ylim(0,5)
#subplot(6,1,6)
#plot(transpose(gi_conductance_over_time[Ne*number_patterns + Ni + 1,:]),c="b")
#ylabel("Inhibitory conductance")
#title("Pattern B, inh neuron 1")
#ylim(-66,-10)

tight_layout()
savefig(string("figures/",datestring,"_conductance_for_trigger_with_two_pat_detector.png"),dpi=300)
=#
