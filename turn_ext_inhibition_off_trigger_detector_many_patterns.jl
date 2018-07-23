# This file is made by Marcelo Gauy and is based on
# simulation code from Litwin-Kumar 2014.

# This simulation checks whether we can turn on and off the external inhibition and get the two pattern
# detector to shut down activity in E[i] 

# Package for plotting
using PyPlot
# This is just a package to save/load from file
using HDF5

using JLD

include("simulation5.jl")
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

number_patterns = 6
#number_patterns = 2
#number_patterns = 3
#number_patterns = 10

function weightpars_custom() #parameters needed to generate weight matrix
	return Ne,Ni,jee0,jei0,jie,jii,jba,jca,jac,p,rex,rix,jex,jix,j_back
end

function stimpars_custom()
    return [1.0 24.0 0.25 30.0]
end

_,_,stim_rate = stimpars_custom()

T = 2000 #6 patterns
#T = 1200 #2 patterns
#T = 6000 #10 patterns
#T = 800
dt,_,Nskip,vpeak,Nspikes = simpars()

function simpars_custom()
    return dt, T, Nskip, vpeak, Nspikes
end


turn_off_times = [200 550 1000 1300 1700] #6 patterns
#turn_off_times = [550] #2 patterns
#turn_off_times = [550 1200]#  3 patterns
time_off = 100
#turn_on_times = T#turn_off_times + time_off
#turn_on_times[3] = T
external_inhibition_turn_off_stim = [turn_off_times[1] T 1.0 0.0; turn_off_times[2] T 1.0 0.0; turn_off_times[3] T 1.0 0.0; turn_off_times[4] T 1.0 0.0; turn_off_times[5] T 1.0 0.0] #6patterns
#external_inhibition_turn_off_stim = [turn_off_times[1] T 1.0 0.0] #2 patterns
#external_inhibition_turn_off_stim = [turn_off_times[1] T 1.0 0.0; turn_off_times[2] T 1.0 0.0] #3patterns
#turn_off_times = [200 550 1000 1300 1700 2000 2500 2900 5600] #10 patterns
#external_inhibition_turn_off_stim = [turn_off_times[1] T 1.0 0.0 ; turn_off_times[2] T 1.0 0.0; turn_off_times[3] T 1.0 0.0; turn_off_times[4] T 1.0 0.0; turn_off_times[5] T 1.0 0.0;  turn_off_times[6] T 1.0 0.0;  turn_off_times[7] T 1.0 0.0;  turn_off_times[8] T 1.0 0.0;  turn_off_times[9] T 1.0 0.0] #10patterns

rextin = 0.05 # khz
jextin = 2.3 + 0.2*noise[6] # ???
function inhibition_1_pars()
	return rextin,jextin,external_inhibition_turn_off_stim#[turn_off_time T 1.0 0.0]
end


println("jee0 =", jee0)
println("jei0 =", jei0)
println("jie =", jie)
println("jii =", jii)
println("jba =", jba)
println("jca =", jca)
println("jac =", jac)
println("jextin =", jextin)


tic()

times,ns,times_x_axis,voltage_over_time,geaff_conductance_over_time,gerec_conductance_over_time,gi_conductance_over_time,spike_count_over_time,turn_on_time, ex_spike_count_per_step_per_pattern, in_spike_count_per_step_per_pattern, times_global_inhibitory_cells, ns_global_inhibitory_cells = sim(weightpars_custom,membranepars,synapsepars,stimpars_custom,simpars_custom,inhibition_1_pars, number_patterns)
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
#smoothing_time = 10 #divisor of T
#rate_over_time_smoothed = zeros(number_patterns, T)
#for pattern_index = 1:number_patterns
#	for time_index = 1:T
#		if (time_index < smoothing_time)
#			aux_spike_count = sum(spike_count_per_step_per_pattern[pattern_index, 1:10*time_index]) # the 10 is 1/dt of the simulation
#			rate_over_time_smoothed[pattern_index,time_index] = 1000*aux_spike_count/(Ne*time_index)
#		else
#			aux_spike_count = sum(spike_count_per_step_per_pattern[pattern_index, (10*time_index-10*smoothing_time+1):10*time_index])
#			rate_over_time_smoothed[pattern_index,time_index] = 1000*aux_spike_count/(Ne*smoothing_time)
#		end
#	end
#end
###################################


#####smoothing out the spike count
smoothing_time = 5 #divisor of T
ex_rate_over_time_smoothed = zeros(number_patterns, T)
in_rate_over_time_smoothed = zeros(number_patterns+2, T) ###last position contains the global inhibitory cells.
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
fig = plt[:figure](figsize=(14,7.5))
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
ax[:set_xlabel]("Time(s)", fontsize=30)
ax[:set_ylabel]("Ensemble", fontsize=30)
ax[:yaxis][:labelpad] = 40
axsec = ax[:twinx]()
axsec[:set_ylabel]("Rate(Hz)", fontsize=30)
axsec[:spines]["top"][:set_color]("none")
axsec[:spines]["bottom"][:set_color]("none")
axsec[:spines]["left"][:set_color]("none")
axsec[:spines]["right"][:set_color]("none")
axsec[:tick_params](labelcolor="w", top="off", bottom="off", left="off", right="off")
#figure(figsize=(12,2*7))
#title("Test in turning off external inhibition with a two pattern detector to prevent both patterns being active.")
#ax = gca()
for pattern_index = 1:number_patterns
	##excitatory raster plots
	rowcount = 1
	#fig[:add_subplot](2*number_patterns + 1, 1, 2*pattern_index - 1, sharex = ax)
	fig[:add_subplot](2*number_patterns + 2, 1, 2*pattern_index - 1, sharex = ax) #2 patterns figure 2D without inhibitory ensembles
	ax1 = gca()
	setp(ax1[:get_xticklabels](),visible=false) # Disable x tick labels
	setp(ax1[:get_yticklabels](),fontsize=15) # Disable x tick labels
	#ax1[:tick_params](axis="y", which="both", left="off", labelleft="off")
		#ax1[:tick_params](labelsize=20)
		#ax1[:set_xticklabels](["0","0","0.5","1.0","1.5","2.0"])#,"2.5","3.0","3.5"])
	for cc = (pattern_index-1)*Ne+1:pattern_index*Ne
		vals = times[cc,1:ns[cc]]
		y = rowcount*ones(length(vals))
		scatter(vals,y,s=2,c="k",marker="o",linewidths=0)
		#ax[:scatter](vals,y,s=.3,c="k",marker="o",linewidths=0)
		rowcount+=1
	end
	if number_patterns == 2
		#plot([turn_off_times[1]; turn_off_times[1]],[0; T],c="k")
		#plot([turn_on_time[2]; turn_on_time[2]],[0; T],c="k")
	else
		for index = 1:size(turn_off_times)[2]## with turn_off_times having more than 1 element
			plot([turn_off_times[index]; turn_off_times[index]],[0; T],c="k")
			if (turn_on_time[index] != 0)
				plot([turn_on_time[index]; turn_on_time[index]],[0; T],c="b")
			end
		end
		plot([turn_on_time[size(turn_off_times)[2]+1]; turn_on_time[size(turn_off_times)[2]+1]],[0; T],c="b")
	end
	xlim(0,T)
	ylim(0,Ne)
	if pattern_index == 1
		#ylabel("\$E_1\$",fontsize=25)
		ylabel("\$E_1\$",fontsize=20) # figure 2D 2 patterns no inhibition in ensembles
	elseif pattern_index == 2
		#ylabel("\$E_2\$",fontsize=25)
		ylabel("\$E_2\$",fontsize=20) # figure 2D 2 patterns no inhibition in ensembles
	elseif pattern_index == 3
		ylabel("\$E_3\$",fontsize=25)
	elseif pattern_index == 4
		ylabel("\$E_4\$",fontsize=25)
	elseif pattern_index == 5
		ylabel("\$E_5\$",fontsize=25)
	elseif pattern_index == 6
		ylabel("\$E_6\$",fontsize=25)
	end
	ax2 = ax1[:twinx]()
	ax2[:plot](ex_rate_over_time_smoothed[pattern_index,:],c="r")
	#ax2[:set_ylabel]("Ex. Rate", color = "r")
	ax2[:tick_params](colors="r")
	ax2[:tick_params](labelsize = 15)
	ax2[:tick_params](axis="y", which="both", left="off")
	ylim(0,80)
	
	##inhibitory raster plots - not in figure 2D
	fig[:add_subplot](2*number_patterns + 2, 1, 2*pattern_index, sharex = ax)  #2 patterns figure 2D without inhibitory ensembles
	ax2 = gca()
	setp(ax2[:get_xticklabels](),visible=false) # Disable x tick labels
	setp(ax2[:get_yticklabels](),fontsize=15) # Disable x tick labels
	rowcount = 1
	for cc = 1+(pattern_index-1)*Ni + number_patterns*Ne:number_patterns*Ne+pattern_index*Ni
		vals = times[cc,1:ns[cc]]
		y = rowcount*ones(length(vals))
		scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
		#ax[1][:scatter](vals,y,s=.3,c="k",marker="o",linewidths=0)
		rowcount+=1
	end
	if number_patterns == 2
	#	plot([turn_off_times[1]; turn_off_times[1]],[0; T],c="k")
	#	plot([turn_on_time[2]; turn_on_time[2]],[0; T],c="b")
	else
		for index = 1:size(turn_off_times)[2]## with turn_off_times having more than 1 element
			plot([turn_off_times[index]; turn_off_times[index]],[0; T],c="k")
			if (turn_on_time[index] != 0)
				plot([turn_on_time[index]; turn_on_time[index]],[0; T],c="b")
			end
		end
		plot([turn_on_time[size(turn_off_times)[2]+1]; turn_on_time[size(turn_off_times)[2]+1]],[0; T],c="b")
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
	elseif pattern_index == 7
		ylabel("\$I_7\$",fontsize=25)
	elseif pattern_index == 8
		ylabel("\$I_8\$",fontsize=25)
	elseif pattern_index == 9
		ylabel("\$I_9\$",fontsize=25)
	elseif pattern_index == 10
		ylabel("\$I_{10}\$",fontsize=25)
	end
	ax3 = ax2[:twinx]()
	ax3[:plot](in_rate_over_time_smoothed[pattern_index,:],c="b")
	#in_rate_over_time_smoothed[pattern_index,:]
	#ax3[:set_ylabel]("In. Rate(Hz)", color = "b")
	ax3[:tick_params](colors="b")
	ax3[:tick_params](labelsize = 15)
	ylim(0,130)
	#xlabel("Time")
end
#pattern detector inhibition raster plots

#fig[:add_subplot](2*number_patterns + 1, 1, 2*number_patterns + 1, sharex = ax)
fig[:add_subplot](2*number_patterns + 2, 1, 2*number_patterns + 1, sharex = ax) # 2 patterns figure 2D without inhibitory ensembles
ax2 = gca()
setp(ax2[:get_xticklabels](),visible=false) # Disable x tick labels
setp(ax2[:get_yticklabels](),fontsize=15) # Disable x tick labels
#ax2[:tick_params](axis="y", which="both", left="off", labelleft="off")
#ax2[:tick_params](labelsize=20)
#ax2[:set_xticklabels](["0","0","0.25","0.5","0.75","1.0","1.25","1.5","1.75","2.0"])
#ax2[:set_xticklabels](["0","0","0.2","0.4","0.6","0.8","1.0","1.5","1.75","2.0"])  # 2 patterns figure 2D without inhibitory ensembles
#setp(ax2[:get_yticklabels](),visible=false) # Disable x tick labels
rowcount = 1
for cc = 1+number_patterns*Ne+number_patterns*Ni:number_patterns*Ne+(number_patterns+1)*Ni
	vals = times[cc,1:ns[cc]]
	y = rowcount*ones(length(vals))
	scatter(vals,y,s=2,c="b",marker="o",linewidths=0)
	rowcount+=1
end
if number_patterns == 2
	#plot([turn_off_times[1]; turn_off_times[1]],[0; T],c="k")
	#plot([turn_on_time[2]; turn_on_time[2]],[0; T],c="k")
else
	for index = 1:size(turn_off_times)[2]## with turn_off_times having more than 1 element
		plot([turn_off_times[index]; turn_off_times[index]],[0; T],c="k")
		if (turn_on_time[index] != 0)
			plot([turn_on_time[index]; turn_on_time[index]],[0; T],c="b")
		end
	end
	plot([turn_on_time[size(turn_off_times)[2]+1]; turn_on_time[size(turn_off_times)[2]+1]],[0; T],c="b")
end
ylabel("\$I_C\$",fontsize=20)
#xlim(0,T)
ylim(0,Ni)
ax3 = ax2[:twinx]()
ax3[:plot](in_rate_over_time_smoothed[number_patterns+1,:],c="b")
#in_rate_over_time_smoothed[number_patterns+1,:]
#ax3[:set_ylabel]("In. Rate", color = "b")
ax3[:tick_params](colors="b")
ax3[:tick_params](labelsize = 15)
#ax3[:tick_params](axis="y", which="both", left="off")
#xlabel("Time")


#fig[:add_subplot](2*number_patterns + 1, 1, 2*number_patterns + 1, sharex = ax)
fig[:add_subplot](2*number_patterns + 2, 1, 2*number_patterns + 2, sharex = ax) # 2 patterns figure 2D without inhibitory ensembles
ax4 = gca()
ax4[:tick_params](labelsize=20)
setp(ax4[:get_yticklabels](),fontsize=15) # Disable x tick labels
#ax2[:tick_params](axis="y", which="both", left="off", labelleft="off")
#ax4[:set_xticklabels](["0","0","0.25","0.5","0.75","1.0","1.25","1.5","1.75","2.0"])
#ax4[:set_xticklabels](["0","0","0.2","0.4","0.6","0.8","1.0","1.5","1.75","2.0"])  # 2 patterns figure 2D without inhibitory ensembles
#setp(ax4[:get_yticklabels](),visible=false) # Disable x tick labels
#rowcount = 1
#for cc = 1+number_patterns*Ne+number_patterns*Ni:number_patterns*Ne+(number_patterns+1)*Ni
#	vals = times[cc,1:ns[cc]]
#	y = rowcount*ones(length(vals))
#	scatter(vals,y,s=2,c="b",marker="o",linewidths=0)
#	rowcount+=1
#end
rowcount = 1
for cc = 1:100
	vals = times_global_inhibitory_cells[cc,1:ns_global_inhibitory_cells[cc]]
	y = rowcount*ones(length(vals))
	scatter(vals,y,s=2,c="b",marker="o",linewidths=0)
	rowcount+=1
end
if number_patterns == 2
#	plot([turn_off_times[1]; turn_off_times[1]],[0; T],c="k")
#	plot([turn_on_time[2]; turn_on_time[2]],[0; T],c="k")
else
	for index = 1:size(turn_off_times)[2]## with turn_off_times having more than 1 element
		plot([turn_off_times[index]; turn_off_times[index]],[0; T],c="k")
		if (turn_on_time[index] != 0)
			plot([turn_on_time[index]; turn_on_time[index]],[0; T],c="b")
		end
	end
	plot([turn_on_time[size(turn_off_times)[2]+1]; turn_on_time[size(turn_off_times)[2]+1]],[0; T],c="b")
end
ylabel("\$I_G\$",fontsize=20)
ax4[:yaxis][:labelpad] = 40
xlim(0,T)
ylim(0,Ni)
ax5 = ax4[:twinx]()
ax5[:plot](in_rate_over_time_smoothed[number_patterns+2,:],c="b")
#ax3[:set_ylabel]("In. Rate", color = "b")
ax5[:tick_params](colors="b")
ax5[:tick_params](labelsize = 15)
ax5[:tick_params](axis="y", which="both", left="off")


tight_layout()
subplots_adjust(hspace=0.0)
savefig(string("figures/",datestring,"spiking_two_patterns_dynamics.png"),dpi=300) #figure 2D

toc()




save("data/replay_data.jld", "times", times, "ns", ns, "times_x_axis", times_x_axis, "turn_on_time", turn_on_time, "turn_off_times", turn_off_times, "ex_spike_count_per_step_per_pattern", ex_spike_count_per_step_per_pattern, "in_spike_count_per_step_per_pattern", in_spike_count_per_step_per_pattern, "times_global_inhibitory_cells", times_global_inhibitory_cells, "ns_global_inhibitory_cells", ns_global_inhibitory_cells)


#=
rowcount = 1
for pattern_index = 1:number_patterns
	for cc = (pattern_index-1)*Ne+1:pattern_index*Ne
		subplot(2*number_patterns+2,1,pattern_index)
		vals = times[cc,1:ns[cc]]
		y = rowcount*ones(length(vals))
		scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
		rowcount+=1
	end
	for index = 1:size(turn_off_times)[2]
		plot([turn_off_times[index]; turn_off_times[index]],[0; T],c="r")
		plot([turn_on_time[index]; turn_on_time[index]],[0; T],c="b")
	end
	plot([turn_on_time[size(turn_off_times)[2]+1]; turn_on_time[size(turn_off_times)[2]+1]],[0; T],c="b")
	ylabel("PE")
	xlim(0,T)
	ylim(0,Ne)
	rowcount = 1
end
for pattern_index = 1:number_patterns
	for cc = 1+(pattern_index-1)*Ni + number_patterns*Ne:number_patterns*Ne+pattern_index*Ni
		subplot(2*number_patterns+2,1,number_patterns + pattern_index)
		vals = times[cc,1:ns[cc]]
		y = rowcount*ones(length(vals))
		scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
		rowcount+=1
	end
	for index = 1:size(turn_off_times)[2]
		plot([turn_off_times[index]; turn_off_times[index]],[0; T],c="r")
		plot([turn_on_time[index]; turn_on_time[index]],[0; T],c="b")
	end
	plot([turn_on_time[size(turn_off_times)[2]+1]; turn_on_time[size(turn_off_times)[2]+1]],[0; T],c="b")
	ylabel("PI")
	xlim(0,T)
	ylim(0,Ni)
	rowcount = 1
end

for cc = 1+number_patterns*Ne+number_patterns*Ni:number_patterns*Ne+(number_patterns+1)*Ni
	subplot(2*number_patterns+2,1,2*number_patterns+1)
	vals = times[cc,1:ns[cc]]
	y = rowcount*ones(length(vals))
	scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
	rowcount+=1
end
for index = 1:size(turn_off_times)[2]
	plot([turn_off_times[index]; turn_off_times[index]],[0; T],c="r")
	plot([turn_on_time[index]; turn_on_time[index]],[0; T],c="b")
end
plot([turn_on_time[size(turn_off_times)[2]+1]; turn_on_time[size(turn_off_times)[2]+1]],[0; T],c="b")
ylabel("PDI")
xlim(0,T)
ylim(0,Ni)

subplot(2*number_patterns+2,1,2*number_patterns+2)
for pattern_index = 1:1
	plot(rate_over_time_smoothed[pattern_index,:],c="b")
end
xlim(0,T)
ylabel("Rate")
xlabel("Time")
#for index = 1:size(turn_off_times)[2]
#	plot([10*turn_off_times[index]; 10*turn_off_times[index]],[0; 4],c="r")
#	plot([10*turn_on_time[index]; 10*turn_on_time[index]],[0; 4],c="b")
#end
tight_layout()
savefig(string("figures/",datestring,"_raster_for_trigger_with_many_pat_detector.png"),dpi=300)

toc()
=#



#clf()
#figure(figsize=(12,2*4))
#subplot(4,1,1)
#plot(transpose(geaff_conductance_over_time[1,:]),c="b")
#plot(transpose(gi_conductance_over_time[1,:]),c="k")
#plot(transpose(geaff_conductance_over_time[Ne+1,:]),c="r")
#ylabel("Excitatory conductance")
#title("Pattern A, neuron 1")
#ylim(0,5)

#subplot(4,1,2)
#plot(transpose(geaff_conductance_over_time[Ne+1,:]),c="b")
#ylabel("Excitatory conductance")
#title("Pattern B, neuron 1")
#ylim(0,5)

#subplot(4,1,3)
#plot(transpose(mean(geaff_conductance_over_time[1:Ne,:],1)),c="b")
#plot(transpose(mean(gi_conductance_over_time[1:Ne,:],1)),c="k")
#plot(transpose(mean(geaff_conductance_over_time[Ne+1:2*Ne,:],1)),c="r")
#plot(transpose(mean(gi_conductance_over_time[Ne+1:2*Ne,:],1)),c="g")
#ylabel("Excitatory conductance")
#title("Pattern A, average")
#ylim(0,5)

#subplot(4,1,4)
#plot(transpose(mean(geaff_conductance_over_time[Ne+1:2*Ne,:],1)),c="b")
#plot(transpose(mean(geaff_conductance_over_time[2*Ne+1:3*Ne,:],1)),c="r")
#plot(transpose(mean(geaff_conductance_over_time[3*Ne+1:4*Ne,:],1)),c="g")
#plot(transpose(mean(geaff_conductance_over_time[4*Ne+1:5*Ne,:],1)),c="k")
	#for index = 1:size(turn_off_times)[2]
	#	plot([10*turn_off_times[index]; 10*turn_off_times[index]],[0; 10],c="r")
	#	plot([10*turn_on_times[index]; 10*turn_on_times[index]],[0; 10],c="b")
	#end
#ylabel("Excitatory conductance")
#title("Pattern B, average")

#tight_layout()
#savefig(string("figures/",datestring,"_conductance_for_trigger_with_many_pat_detector.png"),dpi=300)


#clf()
#figure(figsize=(12,2*4))
#subplot(4,1,1)
#plot(times_x_axis,voltage_over_time[1,:],c="b")
#ylabel("Membrane voltage")
#title("Pattern A, neuron 1")
#ylim(-66,-10)

#subplot(4,1,2)
#plot(times_x_axis,voltage_over_time[2,:],c="b")
#ylabel("Membrane voltage")
#title("Pattern A, neuron 2")
#ylim(-66,-10)

#subplot(4,1,3)
#plot(times_x_axis,voltage_over_time[3,:],c="b")
#ylabel("Membrane voltage")
#title("Pattern B, neuron 1")
#ylim(-66,-10)

#subplot(4,1,4)
#plot(times_x_axis,voltage_over_time[4,:],c="b")
#ylabel("Membrane voltage")
#title("Pattern B, neuron 2")
#ylim(-66,-10)

#tight_layout()
#savefig(string("figures/",datestring,"_voltage_for_trigger_with_two_pat_detector.png"),dpi=300)
