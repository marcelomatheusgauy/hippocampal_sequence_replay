# This file is made by Marcelo Gauy and is based on
# simulation code from Litwin-Kumar 2014.

using Distributions

# 1. We model l balanced units (simulation.jl corresponds to one balanced unit/simulation4.jl corresponds to 2 balanced units)
#    and we refer to them as pattern E[i], 1 <= i <= l.
# 2. E[i] sends excitatory connections to E[i+1] of minimum possible weight (see a_to_b_min_weight.jl to binary search for that weight)
# 3. We send a common external inhibitory signal (parameter inhibition_1_pars) to both a_ex and b_ex such that activation of E[i] should not trigger E[i+1] despite
#    E[i] sending input to E[i+1]. Also, the external inhibition should not kill pattern E[i] as well when it is active.
# 4. If we turn off the external inhibition then E[i+1] can turn active (ext_inhibition_1_min_weight.jl should confirm this)
# 5. In this simulation we have added a two pattern detector (new inhibitory population) which is only activated if both E[i] and E[i+1] are co-active.
# 6. this simulation should have the start be synchronous and activity then proceeds synchronously

include("afferent_recurrent_neuron.jl")

# The file which should be used to call this one is synfiremode.jl

function sim(weightpars,membranepars,synapsepars,stimpars,simpars,inhibition_1_pars, number_patterns) #runs simulation given weight matrix and populations
	###########################################################################
	# Setting parameters
	###########################################################################
	# Network weights
	Ne,Ni,jee0,jei0,jie,jii,jba,jca,jac,p,rex,rix,jex,jix,j_back,ach_factor = weightpars()   
	default_jex = jex
	#number_patterns = 3
	Ncells = number_patterns*Ne + (number_patterns+1)*Ni

	# Setting up the weight matrix
	weights = zeros(Ncells,Ncells)
	# Populations (for clarity of code below)
	c_in = (number_patterns*Ne+number_patterns*Ni+1):Ncells
	for pattern_index = 1:number_patterns
		a_ex = ((pattern_index-1)*Ne+1):(pattern_index*Ne)
		a_in = (number_patterns*Ne+(pattern_index-1)*Ni+1):(number_patterns*Ne+pattern_index*Ni)
		weights[a_ex,a_ex] = jee0*ach_factor
		weights[a_ex,a_in] = jie*ach_factor
		weights[a_in,a_ex] = jei0*1.2
		weights[a_in,a_in] = jii*1.2

		if(pattern_index < number_patterns)
			a_ex_next = ((pattern_index)*Ne+1):((pattern_index+1)*Ne)
			a_in_next = (number_patterns*Ne+(pattern_index)*Ni+1):(number_patterns*Ne+(pattern_index+1)*Ni)
			weights[a_ex,a_ex_next] = jba*ach_factor
			weights[a_ex,a_in_next] = jba*(jie/jee0)*ach_factor
			weights[a_ex_next,a_ex] = jba*j_back*ach_factor
			weights[a_ex_next,a_in] = jba*(jie/jee0)*j_back*ach_factor
		end
		weights[a_ex,c_in] = jca*ach_factor
		weights[c_in,a_ex] = jac*1.2
		weights[c_in,c_in] = jii*1.2
	end

##################################################
	#break the symmetry in the system
#	random_connectivity = rand(Ncells,Ncells) .< p

#	for pattern_index = 2:number_patterns
#		high_density_rand_subset = zeros(Ncells)
#		previous_pattern = zeros(Ncells)
#		previous_pattern[(pattern_index-2)*Ne+1:(pattern_index-1)*Ne] = ones(Ne)
#		high_density_rand_subset[(pattern_index-1)*Ne+1:(pattern_index)*Ne] = (rand(Ne) .< 0.2)
#		extra_connectivity = previous_pattern*transpose(high_density_rand_subset)
#		extra_connectivity_rand = extra_connectivity.*(rand(Ncells,Ncells).<0.2)
#		random_connectivity = random_connectivity + extra_connectivity_rand
#	end
#	weights = weights.*random_connectivity
##################################################


	weights = weights.*(rand(Ncells,Ncells) .< p)
	for cc = 1:Ncells
        	weights[cc,cc] = 0
	end
	Ne *= number_patterns #<= what was this line for???
	Ni *= (number_patterns+ 1) #<= what was this line for??? replaced 2 to number_patterns and 3 for number_patterns+1

############################################################

	# Membrane potential params
	taue,taui,vleake,vleaki,deltathe,C,erev,irev,vth0,ath,tauth,vre,taurefrac,aw_adapt,bw_adapt,tauw_adapt = membranepars()

	# Synapse specific params
	tauerise,tauedecay,tauirise,tauidecay,nmdaspikethreshold,dendriticA,dendriticB,dendriticC,taudendriticA,taudendriticB,taudendriticC,dendriticlatency, dendriticintegrationwindow, dendriticrefractoryperiod = synapsepars()

	# Columns are: (1) Rate, (2) start time, (3) end time
	stims = stimpars()
	stim_rate = 1.0
	default_stim_rate = stim_rate
	stim_weight = 0.0

	# External common inhibitory input
	rextin,jextin,ext_in_rates = inhibition_1_pars()
	default_rextin = rextin
	default_jextin = jextin

	# simulation
	dt,T,Nskip,vpeak,Nspikes = simpars()

	###########################################################################
	# Initializing variables
	###########################################################################
	times = zeros(Ncells,Nspikes)
	ns = zeros(Int,Ncells)

	forwardInputsEaff = zeros(Ncells) #summed weight of incoming afferent E spikes
	forwardInputsEPrevaff = zeros(Ncells) #as above, for previous timestep
	forwardInputsErec = zeros(Ncells) #summed weight of incoming recurrent E spikes
	forwardInputsEPrevrec = zeros(Ncells) #as above, for previous timestep

	forwardInputsI = zeros(Ncells)
	forwardInputsIPrev = zeros(Ncells)

	xeriseaff = zeros(Ncells) #auxiliary variables for afferent E/I currents (difference of exponentials)
	xedecayaff = zeros(Ncells)
	xeriserec = zeros(Ncells) #auxiliary variables for recurrent E/I currents (difference of exponentials)
	xedecayrec = zeros(Ncells)

	xirise = zeros(Ncells)
	xidecay = zeros(Ncells)

	##############################dendritic spikes
	###stereotypical current pulse of a dendritic spike using physiological parameters
	###we only take the values it assumes in the first 1ms (this is the region where it is positive)
	###in the code there is a 2.7ms latency delay between the moment the dendritic spike is sent forward
	###dendritic spikes are initiated if conductance changes drastically in a 2ms interval
	###there is a 5ms refractory period on dendritic spikes
	exponentialdendriticA = dendriticA*exp(-(0:30)/taudendriticA)
	exponentialdendriticB = dendriticB*exp(-(0:30)/taudendriticB)
	exponentialdendriticC = dendriticC*exp(-(0:30)/taudendriticC)
	dendriticInputcurve = -exponentialdendriticA+exponentialdendriticB-exponentialdendriticC 

	dendriticspiketimes = -50*ones(Ncells)
	dendriticspiketimes[1:200] = 100 #activates the first pattern

	expdist = Exponential()

	v = zeros(Ncells) #membrane voltage 
	nextx = zeros(Ncells) #time of next external excitatory input
	nextin = zeros(Ncells) #time of next external inhibitory input (only onto excitatory neurons)
	nextxstim = zeros(Ncells) # time of next stimulus input
	rx = zeros(Ncells) #rate of external input

	for cc = 1:Ncells
		#v[cc] = vre + (vth0-vre)*rand() #might make sense for the synchronous mode to change this -> this is the reason why the initial pattern does not spike so synchronously now
		v[cc] = vre # to make first pattern synchronous
		# Next external spike
		if cc <= Ne 
			rx[cc] = rex
			nextx[cc] = rand(expdist)/rx[cc]
		else
			rx[cc] = rix
			nextx[cc] = rand(expdist)/rx[cc]
		end
	end

	# External inhibition spike
	for cc = 1:Ncells
		#if cc <= Ncells
		nextin[cc] = rand(expdist)/rextin
		#end
	end

	vth = vth0*ones(Ncells) #adaptive threshold
	wadapt = aw_adapt*(vre-vleake)*ones(Ne) #adaptation current
	lastSpike = -100*ones(Ncells) #last time the neuron spiked

	Nsteps = round(Int,T/dt)
	spike_count_per_step_per_pattern = zeros(number_patterns,Nsteps)
	spike_count_over_time_per_pattern = zeros(number_patterns,Nsteps)
	ex_spike_count_per_step_per_pattern = zeros(number_patterns,Nsteps)
	in_spike_count_per_step_per_pattern = zeros(number_patterns + 1,Nsteps)


	voltage_over_time = zeros(4,Nsteps) # To keep track of voltage
	geaff_conductance_over_time = zeros(Ncells,Nsteps) # To keep track of voltage
	gerec_conductance_over_time = zeros(Ncells,Nsteps) # To keep track of voltage
	gi_conductance_over_time = zeros(Ncells,Nsteps) # To keep track of voltage
	time_xaxis = zeros(1,Nsteps)
	turn_on_time = zeros(number_patterns)
	auxiliary_pattern_index = 1
	inhibitory_stimulus_off = false


	###########################################################################
	# Main simulation loop
	###########################################################################
	for tt = 1:Nsteps
		if mod(tt,Nsteps/100) == 1  #print percent complete
			@printf("\r%d%%",round(Int,100*tt/Nsteps))
		end
		t = dt*tt
		time_xaxis[1,tt] = t

		#voltage_over_time[1,tt] = v[1]
		#voltage_over_time[2,tt] = v[2]
		#voltage_over_time[3,tt] = v[1+round(Int,Ne/2)]
		#voltage_over_time[4,tt] = v[2+round(Int,Ne/2)]
		for pattern_index in 1:number_patterns
			voltage_over_time[1,tt] = v[1+round(Int, (Ne*(pattern_index-1))/number_patterns)]######
			voltage_over_time[2,tt] = v[2+round(Int, (Ne*(pattern_index-1))/number_patterns)]######
		end

		# These variables collect new input spikes for the neuron
		forwardInputsEaff[:] = 0.
		forwardInputsErec[:] = 0.

		forwardInputsI[:] = 0.

		# update single cells
		spiked = zeros(Bool,Ncells)	

		tprev = dt*(tt-1)
		for cc = 1:Ncells

			# Checking if a stimulus is on
            		if cc <= Ne/number_patterns  || Ne+1 <= cc <= Ne + Ni/(number_patterns+1)
              			foundstim = false
				for i in 1:size(stims)[1]
                   			if stims[i,1] <= t <= stims[i,2]
						if cc <= Ne/number_patterns
                            				stim_rate = stims[i,3]
                            				stim_weight = stims[i,4]
                        			else
	                            			stim_rate = stims[i,3]
	                            			stim_weight = stims[i,4] * jie/jee0
                        			end
                   				foundstim = true
                    			end
				end
				# Set to background rate
				if !foundstim
                    			stim_rate = default_stim_rate
                   			stim_weight = 0.0
				end
            		end

			# Updating external inhibition rates
			#if cc <= Ne
			#foundstim = false
			for i in 1:size(ext_in_rates)[1]
				#if ext_in_rates[i,1] <= t <= ext_in_rates[i,2]#########adapted to stop artificially
				if ext_in_rates[i,1] == t
					rextin = ext_in_rates[i,3]
					jextin = ext_in_rates[i,4]
					if rextin <= 0.0001 && nextin[cc] < ext_in_rates[i,2]
						nextin[cc] = ext_in_rates[i,2] + rand(expdist)/default_rextin
						#println(nextin[cc])
					end
					#println(rextin)
					inhibitory_stimulus_off = true #########adapted to stop artificially
				end
			end
			# Set to normal background rate
			#if !foundstim
			#	rextin = default_rextin
			#	jextin = default_jextin
			#end
			#end

			#############################################artificially turns on stopping inhibition -> synfire does not have it
			#if(inhibitory_stimulus_off)
			#	if(tt > 380)
			#		if (spike_count_over_time_per_pattern[auxiliary_pattern_index+1,tt - 380] > 10) #this number might have to be adapted with the system
			#			auxiliary_pattern_index += 1
			#			turn_on_time[auxiliary_pattern_index] = t;
			#			rextin = default_rextin
			#			jextin = default_jextin
			#			inhibitory_stimulus_off = false
			#			println(auxiliary_pattern_index)
			#		end
			#	end
			#end
			#########################################


			# External input spike generation
			while(t > nextx[cc]) 
				nextx[cc] += rand(expdist)/rx[cc]
				if cc <= Ne
                   			if jex > 0
						forwardInputsEPrevrec[cc] += jex
                    			else
						forwardInputsIPrev[cc] -= jex
					end
				else
					if jix > 0
						forwardInputsEPrevrec[cc] += jix
					else
						forwardInputsIPrev[cc] -= jix
					end
				end
			end

			# Stimulus input generation
            		if (1 <= cc <= Ne/number_patterns || Ne + 1 <= cc <= Ne + Ni/(100*number_patterns+1))
				if (t == 5)
					nextxstim[cc] += 1.0
					#println(nextxstim[cc])
	            			if stim_weight > 0
						if (cc <= Ne)
							forwardInputsEPrevrec[cc] += 0 ###unecessary to give stimulus: currently the first pattern is activated through dendritic spikes at 10ms
						else
							forwardInputsEPrevrec[cc] += 0 #inhibitory neurons are activated more easily than excitatory thus the difference
						end
			        	else
	                			forwardInputsIPrev[cc] -= 0 #originally 600
					end
				end
 			end

			# Stimulus input generation #we adapt here for reversal replay
            		if (Ne + Ni - Ni/(number_patterns+1) +1 <= cc <= Ne + Ni) #start at last pattern and not at first...
				if (t == 2)
					nextxstim[cc] += 1.0
					#println(nextxstim[cc])
	            			if stim_weight > 0
						forwardInputsEPrevrec[cc] += 75 #inhibitory neurons are activated more easily than excitatory thus the difference
			        	else
	                			forwardInputsIPrev[cc] -= 0 #originally 600
					end
				end
 			end

			# External inhibitory input onto excitatory neurons
			while(t > nextin[cc])
				if cc <= Ne
					forwardInputsIPrev[cc] += jextin
				else
					forwardInputsIPrev[cc] += jextin * jii/jei0
				end
				nextin[cc] += rand(expdist)/rextin
			end
			
			# Synaptic rise/decay update
			xeriseaff[cc] += -dt*xeriseaff[cc]/tauerise + forwardInputsEPrevaff[cc]
			xedecayaff[cc] += -dt*xedecayaff[cc]/tauedecay + forwardInputsEPrevaff[cc]
			# Synaptic rise/decay update
			xeriserec[cc] += -dt*xeriserec[cc]/tauerise + forwardInputsEPrevrec[cc]
			xedecayrec[cc] += -dt*xedecayrec[cc]/tauedecay + forwardInputsEPrevrec[cc]
			
			xirise[cc] += -dt*xirise[cc]/tauirise + forwardInputsIPrev[cc]
			xidecay[cc] += -dt*xidecay[cc]/tauidecay + forwardInputsIPrev[cc]



			#############HACK################
			#if (lastSpike[cc] <= t <= lastSpike[cc] + taurefrac)
			#	xeriseaff[cc] = 0 
			#	xedecayaff[cc] = 0
			#	xeriserec[cc] = 0 
			#	xedecayrec[cc] = 0
			#	xirise[cc] = 0
			#	xidecay[cc] = 0
			#end
			#################################

			# Adaptation current for excitatory neurons
			if cc <= Ne
				vth[cc] += dt*(vth0 - vth[cc])/tauth;
				wadapt[cc] += dt*(aw_adapt*(v[cc]-vleake) - wadapt[cc])/tauw_adapt;
			end

			geaff = (xedecayaff[cc] - xeriseaff[cc])/(tauedecay - tauerise);

			geaff_conductance_over_time[cc,tt] = geaff;
			gerec = (xedecayrec[cc] - xeriserec[cc])/(tauedecay - tauerise);
			gerec_conductance_over_time[cc,tt] = gerec;
			gi = (xidecay[cc] - xirise[cc])/(tauidecay - tauirise);
			gi_conductance_over_time[cc,tt] = gi;

			# Not in refractory period
			if t > (lastSpike[cc] + taurefrac) 
				# Update excitatory and inhibitory conductance
				geaff = (xedecayaff[cc] - xeriseaff[cc])/(tauedecay - tauerise);

###############################################################################################################################
				#if (cc == 1)
				#	println(dendriticspiketimes[cc])
				#end
				geaff_conductance_over_time[cc,tt] = geaff;
				if (tt > dendriticspiketimes[cc]+dendriticrefractoryperiod+dendriticlatency+30 && tt > dendriticintegrationwindow)
					geaff_slope = geaff_conductance_over_time[cc,tt] - geaff_conductance_over_time[cc,tt-dendriticintegrationwindow]
					if (geaff_slope >= nmdaspikethreshold)
						dendriticspiketimes[cc] = tt
					end
				end
				if (dendriticspiketimes[cc]+dendriticlatency<=tt<dendriticspiketimes[cc]+30+dendriticlatency)#20 is the length of the dendritic spike
					dendritic_time_index = convert(Int64,tt-dendriticspiketimes[cc]-dendriticlatency) + 1
					geaff = 4*(10^3)*dendriticInputcurve[dendritic_time_index]/(erev-v[cc])
					#println(geaff)
				end ###
				if (dendriticspiketimes[cc]>0 && dendriticspiketimes[cc]+dendriticlatency+30<= tt <= dendriticspiketimes[cc]+dendriticrefractoryperiod+dendriticlatency+20)
					geaff = 0
				end
###############################################################################################################################

				gerec = (xedecayrec[cc] - xeriserec[cc])/(tauedecay - tauerise);

				gi = (xidecay[cc] - xirise[cc])/(tauidecay - tauirise);

				# Membrane voltage update
				if cc <= Ne  # Excitatory neuron (eif), has adaptation
					dv = (vleake - v[cc] + deltathe*exp((v[cc]-vth[cc])/deltathe))/taue + geaff*(erev-v[cc])/C + gerec*(erev-v[cc])/C + gi*(irev-v[cc])/C - wadapt[cc]/C;
					v[cc] += dt*dv;
					if v[cc] > vpeak
						spiked[cc] = true
					end
				else
					dv = (vleaki - v[cc])/taui + geaff*(erev-v[cc])/C + gerec*(erev-v[cc])/C + gi*(irev-v[cc])/C;
					v[cc] += dt*dv;
					if v[cc] > vth0
						spiked[cc] = true
					end
				end

				# Spike occurred
				if spiked[cc]
					v[cc] = vre;
					lastSpike[cc] = t;
					ns[cc] += 1;
					#if cc <= Ne
					#	wadapt[cc] += bw_adapt;
					#end
					#xeriseaff[cc] = 0 
					#xedecayaff[cc] = 0
					#xeriserec[cc] = 0 
					#xedecayrec[cc] = 0
					#xirise[cc] = 0
					#xidecay[cc] = 0
					if ns[cc] <= Nspikes
						times[cc,ns[cc]] = t;
					end

					#########################################counting spikes per step to measure rate
					for pattern_index = 1:(number_patterns + 1)
						if (pattern_index < number_patterns + 1)
							if((pattern_index-1)*Ne/number_patterns + 1 <= cc <= pattern_index*Ne/number_patterns)
								ex_spike_count_per_step_per_pattern[pattern_index,tt] += 1
							end
						end
						if(Ne+(pattern_index-1)*Ni/(number_patterns+1) + 1 <= cc <= Ne+pattern_index*Ni/(number_patterns+1))
							in_spike_count_per_step_per_pattern[pattern_index,tt] += 1
						end
					end
					for pattern_index = 1:number_patterns
						if((pattern_index-1)*Ne/number_patterns + 1 <= cc <= pattern_index*Ne/number_patterns)
							spike_count_per_step_per_pattern[pattern_index,tt] += 1
						end
						if(tt == 50)
							for time_index = 1:50
								spike_count_over_time_per_pattern[pattern_index,tt] += spike_count_per_step_per_pattern[pattern_index,tt - time_index + 1]
							end
						elseif (tt > 50)
							spike_count_over_time_per_pattern[pattern_index,tt] = spike_count_over_time_per_pattern[pattern_index,tt-1] + spike_count_per_step_per_pattern[pattern_index,tt] - spike_count_per_step_per_pattern[pattern_index,tt-50]
						end
					end
					#############################################

					if cc <= Ne
						vth[cc] = vth0 + ath;
					end
					
					# Loop over synaptic projections 
					for dd = 1:Ncells
						if cc <= Ne # Excitatory synapse
							if (is_aff(cc,dd,Ne/number_patterns,number_patterns))
								forwardInputsEaff[dd] += weights[cc,dd]; 
							else
								forwardInputsErec[dd] += weights[cc,dd];
							end
						else # Inhibitory synapse
							forwardInputsI[dd] += weights[cc,dd];
						end
					end
				end #end if(spiked)
			end #end if(not refractory)
		end #end loop over cells
		forwardInputsEPrevaff = copy(forwardInputsEaff)
		forwardInputsEPrevrec = copy(forwardInputsErec)
		forwardInputsIPrev = copy(forwardInputsI)
	end #end loop over time
	@printf("\r")

	#println("\n -----------\n")
	#println(size(times))
	#println(maximum(ns))
	if maximum(ns) > Nspikes
		println("Warning! Maximum number of spikes observed was ",maximum(ns)," whereas memory was only allocated for ",Nspikes," spikes.")
	end
	times = times[:,1:maximum(ns)]

	return times,ns,time_xaxis,voltage_over_time,geaff_conductance_over_time,gerec_conductance_over_time,gi_conductance_over_time,spike_count_over_time_per_pattern,turn_on_time, ex_spike_count_per_step_per_pattern,in_spike_count_per_step_per_pattern
end
