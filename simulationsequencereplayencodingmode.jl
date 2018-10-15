# This file is made by Marcelo Gauy and is based on
# simulation code from Litwin-Kumar 2014.

using Distributions

# 1. We model l balanced units (simulation.jl corresponds to one balanced unit/simulation4.jl corresponds to 2 balanced units)
#    and we refer to them as pattern E[i], 1 <= i <= l.
# 2. E[i] sends excitatory connections to E[i+1] of minimum possible weight (see a_to_b_min_weight.jl to binary search for that weight)
# 3. We send a common external inhibitory signal (parameter inhibition_1_pars) to both a_ex and b_ex such that activation of E[i] should not trigger E[i+1] despite
#    E[i] sending input to E[i+1]. Also, the external inhibition should not kill pattern E[i] as well when it is active.
# 4. If we turn off the external inhibition then E[i+1] can turn active
# 5. In this simulation we have added a two pattern detector (new inhibitory population) which is only activated if both E[i] and E[i+1] are co-active.
# 6. In this simulation, we have added place cells; this is a 'simple step', but it is necessary to verify whether sequence ensembles and place ensembles can be associated; This code essentially expands on simulation5.jl by adding parts which include the place cells.

include("afferent_recurrent_neuron.jl")


#we fix the random seed
srand(0)


function sim(weightpars,membranepars,synapsepars,stimpars,simpars,inhibition_1_pars, number_patterns, meanspeed,variancespeed, mousestoppingtimes,mousemovingtimes, cellpathsequence, placecells_per_position) #runs simulation given weight matrix and populations
	###########################################################################
	# Setting parameters
	###########################################################################
	# Network weights
	Ne,Ni,jee0,jei0,jie,jii,jba,jca,jac,p,rex,rix,jex,jix,j_back, N_placecells,rex_place,jex_place,initial_weight_sequence_place = weightpars()   
	default_jex = jex
	Ncells = number_patterns*Ne + (number_patterns+1)*Ni


	############################################################################
	# Plasticity parameters
	
	jeemin = 0.00001 
	jeemax = 5.0 

	##############################################################################

	# Setting up the weight matrix
	weights = zeros(Ncells,Ncells)
	# Populations (for clarity of code below)
	c_in = (number_patterns*Ne+number_patterns*Ni+1):Ncells
	for pattern_index = 1:number_patterns
		a_ex = ((pattern_index-1)*Ne+1):(pattern_index*Ne)
		a_in = (number_patterns*Ne+(pattern_index-1)*Ni+1):(number_patterns*Ne+pattern_index*Ni)
		weights[a_ex,a_ex] = jee0
		weights[a_ex,a_in] = jie
		weights[a_in,a_ex] = jei0
		weights[a_in,a_in] = jii

		if(pattern_index < number_patterns)
			a_ex_next = ((pattern_index)*Ne+1):((pattern_index+1)*Ne)
			a_in_next = (number_patterns*Ne+(pattern_index)*Ni+1):(number_patterns*Ne+(pattern_index+1)*Ni)
			weights[a_ex,a_ex_next] = jba
			weights[a_ex,a_in_next] = jba*(jie/jee0)
			weights[a_ex_next,a_ex] = jba*j_back
			weights[a_ex_next,a_in] = jba*(jie/jee0)*j_back
		end
		weights[a_ex,c_in] = jca
		weights[c_in,a_ex] = jac
		weights[c_in,c_in] = jii

	end


	weights = weights.*(rand(Ncells,Ncells) .< p)

	weights = weights.*(rand(Normal(1, 0.2), Ncells, Ncells))

	for cc = 1:Ncells
        	weights[cc,cc] = 0
	end

	###########################################################################
	# weights to place
	weights_to_place = initial_weight_sequence_place*ones(Ne*number_patterns, N_placecells)

	for index_pat = 1:number_patterns
		pattern_values = (index_pat-1)*200+1:index_pat*200
		weights_to_place[pattern_values[randperm(200)[1:160]],:] = 0
	end

	Ne *= number_patterns
	Ni *= (number_patterns+ 1) 

	############################################################

	# Membrane potential params
	taue,taui,vleake,vleaki,deltathe,C,erev,irev,vth0,ath,tauth,vre,taurefrac,aw_adapt,bw_adapt,tauw_adapt = membranepars()

	# Synapse specific params
	tauerise,tauedecay,tauirise,tauidecay,nmdaspikethreshold,dendriticA,dendriticB,dendriticC,taudendriticA,taudendriticB,taudendriticC,dendriticlatency, dendriticintegrationwindow, dendriticrefractoryperiod = synapsepars()

	# Plasticity parameters
	altd,altp,thetaltd,thetaltp,tauu,tauv,taux = voltagestdppars()

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
	# Initializing variables - sequence
	###########################################################################
	times = zeros(Ncells,Nspikes)
	ns = zeros(Int,Ncells)


	global_inhibitory_cells = 200
	times_global_inhibitory_cells = zeros(global_inhibitory_cells, Nspikes)
	ns_global_inhibitory_cells = zeros(Int, global_inhibitory_cells)

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

	expdist = Exponential()

	v = zeros(Ncells) #membrane voltage 
	nextx = zeros(Ncells) #time of next external excitatory input
	nextin = zeros(Ncells) #time of next external inhibitory input (only onto excitatory neurons)
	nextxstim = zeros(Ncells) # time of next stimulus input
	rx = zeros(Ncells) #rate of external input

	for cc = 1:Ncells
		v[cc] = vre + (vth0-vre)*rand()
		# Next external spike
		if cc <= Ne 
			rx[cc] = rex
			nextx[cc] = rand(expdist)/rx[cc]
		else
			rx[cc] = rix
			nextx[cc] = rand(expdist)/rx[cc]
		end
	end


	number_targets_per_global_inhibitory_cell = zeros(Int, global_inhibitory_cells)
	number_global_inhibitory_targets_per_neuron = floor(Int, global_inhibitory_cells*p/2)
	global_inhibitory_excitatory_inhibitory_targets_to_cell = zeros(Int, Ncells, global_inhibitory_cells)
	
	for cc = 1:Ncells
		auxiliary_targets = randperm(global_inhibitory_cells)[1:number_global_inhibitory_targets_per_neuron]
		global_inhibitory_excitatory_inhibitory_targets_to_cell[cc, auxiliary_targets] = 1
		number_targets_per_global_inhibitory_cell[auxiliary_targets] += 1
	end

	max_number_targets_per_global_inhibitory_cell = maximum(number_targets_per_global_inhibitory_cell)

	global_inhibitory_excitatory_targets = []
	global_inhibitory_inhibitory_targets = []

	auxiliary_excitatory_number_targets_per_global_inhibitory_cell = zeros(Int, global_inhibitory_cells)
	auxiliary_inhibitory_number_targets_per_global_inhibitory_cell = zeros(Int, global_inhibitory_cells)

	global_inhibitory_excitatory_target_per_global_inhibitory_cell = zeros(Int, max_number_targets_per_global_inhibitory_cell)
	global_inhibitory_inhibitory_target_per_global_inhibitory_cell = zeros(Int, max_number_targets_per_global_inhibitory_cell)


	for dd = 1:global_inhibitory_cells
		for cc = 1:Ncells
			if cc < Ne
				if global_inhibitory_excitatory_inhibitory_targets_to_cell[cc,dd] == 1
					auxiliary_excitatory_number_targets_per_global_inhibitory_cell[dd] += 1
					global_inhibitory_excitatory_target_per_global_inhibitory_cell[auxiliary_excitatory_number_targets_per_global_inhibitory_cell[dd]] = cc
				end
			else
				if global_inhibitory_excitatory_inhibitory_targets_to_cell[cc,dd] == 1
					auxiliary_inhibitory_number_targets_per_global_inhibitory_cell[dd] += 1
					global_inhibitory_inhibitory_target_per_global_inhibitory_cell[auxiliary_inhibitory_number_targets_per_global_inhibitory_cell[dd]] = cc
				end
			end
		end
		global_inhibitory_excitatory_target_per_global_inhibitory_cell_aux=global_inhibitory_excitatory_target_per_global_inhibitory_cell[1:auxiliary_excitatory_number_targets_per_global_inhibitory_cell[dd]]
		global_inhibitory_inhibitory_target_per_global_inhibitory_cell_aux=global_inhibitory_inhibitory_target_per_global_inhibitory_cell[1:auxiliary_inhibitory_number_targets_per_global_inhibitory_cell[dd]]

		push!(global_inhibitory_excitatory_targets,global_inhibitory_excitatory_target_per_global_inhibitory_cell_aux)
		push!(global_inhibitory_inhibitory_targets,global_inhibitory_inhibitory_target_per_global_inhibitory_cell_aux)
	end


	vth = vth0*ones(Ncells) #adaptive threshold
	wadapt = aw_adapt*(vre-vleake)*ones(Ne) #adaptation current
	lastSpike = -100*ones(Ncells) #last time the neuron spiked
	u_vstdp = vre*zeros(Ne)
	v_vstdp = vre*zeros(Ne)
	x_vstdp = zeros(Ne)

	Nsteps = round(Int,T/dt)
	dtnormalize = 10 #how often to normalize rows of ee weights
	inormalize = round(Int,dtnormalize/dt)

	ex_spike_count_per_step_per_pattern = zeros(number_patterns,Nsteps)
	in_spike_count_per_step_per_pattern = zeros(number_patterns + 1,Nsteps)
	spike_count_per_step_per_pattern = zeros(number_patterns,Nsteps)
	spike_count_over_time_per_pattern = zeros(number_patterns,Nsteps)

	geaff_conductance_for_nmda_spikes = zeros(Ncells,dendriticintegrationwindow+1)
	time_xaxis = zeros(Nsteps)
	turn_on_time = zeros(number_patterns)
	auxiliary_pattern_index = 1
	auxiliary_pattern_index_place = 1
	inhibitory_stimulus_off = false
	

	###########################################################################
	# Initializing variables - place cells
	###########################################################################
	times_place = zeros(N_placecells,Nspikes)
	ns_place = zeros(Int,N_placecells)

	forwardInputsEplace = zeros(N_placecells) #summed weight of incoming place E spikes
	forwardInputsEPrevplace = zeros(N_placecells) #as above, for previous timestep

	xeriseplace = zeros(N_placecells) #auxiliary variables for afferent E/I currents (difference of exponentials)
	xedecayplace = zeros(N_placecells)

	Normaldist = Normal(meanspeed,variancespeed) # this defines how long the mouse stays at the current position - directly related to the speed it takes to go forward.

	v_place = zeros(N_placecells) #membrane voltage 
	nextx_place = zeros(N_placecells) #time of next external excitatory input
	rx_place = zeros(N_placecells) #rate of external input
	

	for cc = 1:N_placecells
		v_place[cc] = vre + (vth0-vre)*rand()
	end
	
	for index = 1:size(cellpathsequence[:,1],1)
		cc = cellpathsequence[index,1]
		# Next external spike
		rx_place[cc] = rex_place
		nextx_place[cc] = rand(expdist)/rex_place
	end

	vth_place = vth0*ones(N_placecells) #adaptive threshold
	wadapt_place = aw_adapt*(vre-vleake)*ones(N_placecells) #adaptation current
	lastSpike_place = -100*ones(N_placecells) #last time the neuron spiked
	u_vstdp_place = vre*zeros(N_placecells)
	v_vstdp_place = vre*zeros(N_placecells)
	x_vstdp_place = zeros(N_placecells)



	cellpathsequence_length = size(cellpathsequence,2)
	cellpathsequence_times = zeros(cellpathsequence_length,1)
	auxiliar_nextcellpathsequencetime = rand(Normaldist,1)[1]
	index_cellpathsequence = 1
	cellpathsequence_times[index_cellpathsequence] = auxiliar_nextcellpathsequencetime

	cellpathsequence_index = 1
	mouse_stopping_times_index = 1
	external_inhibition_turn_off_stim_index = 1
	mouse_stopped = false

	turn_off_times = zeros(100,1)

	ext_in_rates[external_inhibition_turn_off_stim_index,1] = floor(cellpathsequence_times[cellpathsequence_index]) + 1
	ext_in_rates[external_inhibition_turn_off_stim_index,2] = T
	ext_in_rates[external_inhibition_turn_off_stim_index,3] = 1.0
	ext_in_rates[external_inhibition_turn_off_stim_index,4] = 0.0
	external_inhibition_turn_off_stim_index += 1
	cellpathsequence_index += 1
	println("extin = ",ext_in_rates[external_inhibition_turn_off_stim_index-1,1])
	turn_off_times[external_inhibition_turn_off_stim_index] = ext_in_rates[external_inhibition_turn_off_stim_index,1]

	auxiliary_pattern_recovery_count = zeros(number_patterns)
	auxiliary_time = zeros(number_patterns)

	###########################################################################
	# Main simulation loop
	###########################################################################
	for tt = 1:Nsteps
		if mod(tt,Nsteps/100) == 1  #print percent complete
			@printf("\r%d%%",round(Int,100*tt/Nsteps))
		end
		t = dt*tt
		time_xaxis[tt] = t


		###################################################################
		# go over place cells
		###################################################################
		# These variables collect new input spikes for the neuron
		forwardInputsEplace[:] = 0.


		# update single cells
		spiked_place = zeros(Bool,N_placecells)	

		tprev = dt*(tt-1)
		
		if tt > 10*auxiliar_nextcellpathsequencetime
			auxiliar_nextcellpathsequencetime += rand(Normaldist,1)[1]
			index_cellpathsequence +=1
			cellpathsequence_times[index_cellpathsequence] = auxiliar_nextcellpathsequencetime
			println("auxiliar_nextcellpathsequencetime = ",auxiliar_nextcellpathsequencetime)
			println("index_cellpathsequence = ",index_cellpathsequence)
			println("cellpathsequence_index = ",cellpathsequence_index)
			println("auxiliary_pattern_index = ",auxiliary_pattern_index)
			println("auxiliary_pattern_index_place = ",auxiliary_pattern_index_place)
			println("mousestops =", mousestoppingtimes[mouse_stopping_times_index])
			println("mousemoves =", mousemovingtimes[mouse_stopping_times_index])
			#println("mouse_stopped =", mouse_stopped)

			if (mouse_stopped)
				if(mousemovingtimes[mouse_stopping_times_index] == index_cellpathsequence)
					mouse_stopped = false
					mouse_stopping_times_index += 1
					cellpathsequence_index = index_cellpathsequence
					#println("mousemoves =", mousemovingtimes[mouse_stopping_times_index-1])
					ext_in_rates[external_inhibition_turn_off_stim_index,1] = floor(cellpathsequence_times[cellpathsequence_index])+1
					ext_in_rates[external_inhibition_turn_off_stim_index,2] = T
					ext_in_rates[external_inhibition_turn_off_stim_index,3] = 1.0
					ext_in_rates[external_inhibition_turn_off_stim_index,4] = 0.0
					external_inhibition_turn_off_stim_index += 1
					cellpathsequence_index += 1
				end
			else
				if(mousestoppingtimes[mouse_stopping_times_index] == index_cellpathsequence)
					mouse_stopped = true
					#println("mousestops =", mousestoppingtimes[mouse_stopping_times_index])
					#println("mousemoves =", mousemovingtimes[mouse_stopping_times_index])
					auxiliary_pattern_index_place += 1
				else
					if (index_cellpathsequence == cellpathsequence_index)
						ext_in_rates[external_inhibition_turn_off_stim_index,1] = floor(cellpathsequence_times[cellpathsequence_index])+1
						ext_in_rates[external_inhibition_turn_off_stim_index,2] = T
						ext_in_rates[external_inhibition_turn_off_stim_index,3] = 1.0
						ext_in_rates[external_inhibition_turn_off_stim_index,4] = 0.0
						external_inhibition_turn_off_stim_index += 1
						cellpathsequence_index += 1
						auxiliary_pattern_index_place += 1
						#println("extin = ",ext_in_rates[external_inhibition_turn_off_stim_index-1,1])
						turn_off_times[external_inhibition_turn_off_stim_index] = ext_in_rates[external_inhibition_turn_off_stim_index,1]
					end
				end
			end
			###########################################################################


		end

		if (auxiliary_pattern_index_place == auxiliary_pattern_index)
			for index = 1:size(cellpathsequence[:,index_cellpathsequence],1)
				for placecell_pos_index = 1:placecells_per_position 
					cc = convert(Int64,(placecell_pos_index-1)*N_placecells/placecells_per_position) + cellpathsequence[index,index_cellpathsequence]

					# External input spike generation
					while(t > nextx_place[cc])
						nextx_place[cc] += rand(expdist)/rex_place
						forwardInputsEPrevplace[cc] += jex_place
					end
					# plasticity variables
					u_vstdp_place[cc] += dt*(v_place[cc] - u_vstdp_place[cc])/tauu;
					v_vstdp_place[cc] += dt*(v_place[cc] - v_vstdp_place[cc])/tauv;
					x_vstdp_place[cc] -= dt*x_vstdp_place[cc]/taux;


					# Synaptic rise/decay update
					xeriseplace[cc] += -dt*xeriseplace[cc]/tauerise + forwardInputsEPrevplace[cc]
					xedecayplace[cc] += -dt*xedecayplace[cc]/tauedecay + forwardInputsEPrevplace[cc]

			    		# Adaptation current for excitatory neurons
					vth_place[cc] += dt*(vth0 - vth_place[cc])/tauth;
					wadapt_place[cc] += dt*(aw_adapt*(v_place[cc]-vleake) - wadapt_place[cc])/tauw_adapt;

					geplace = (xedecayplace[cc] - xeriseplace[cc])/(tauedecay - tauerise);

			   	 	# Not in refractory period
					if t > (lastSpike_place[cc] + taurefrac) 
						# Update excitatory and inhibitory conductance
						geplace = (xedecayplace[cc] - xeriseplace[cc])/(tauedecay - tauerise);


			       			# Membrane voltage update
						dv = (vleake - v_place[cc] + deltathe*exp((v_place[cc]-vth_place[cc])/deltathe))/taue + geplace*(erev-v_place[cc])/C - wadapt_place[cc]/C;
						v_place[cc] += dt*dv;
						if v_place[cc] > vpeak
							spiked_place[cc] = true
						end

				 		# Spike occurred
						if spiked_place[cc]
							v_place[cc] = vre;
							lastSpike_place[cc] = t;
							ns_place[cc] += 1;
							if ns_place[cc] <= Nspikes
								times_place[cc,ns_place[cc]] = t;
							end

							if cc <= N_placecells
								vth_place[cc] = vth0 + ath;
							end
						end #end if(spiked)
					end #end if(not refractory)

					##########################################################################################
					# plasticity sequence to place cells
					##########################################################################################

					#vstdp, ltp component
					if (v_place[cc] > thetaltp) && (v_vstdp_place[cc] > thetaltd)
						for dd = 1:Ne
							if weights_to_place[dd,cc] == 0.
								continue
							end

							weights_to_place[dd,cc] += dt*altp*x_vstdp[dd]*(v_place[cc] - thetaltp)*(v_vstdp_place[cc] - thetaltd);
							if weights_to_place[dd,cc] > jeemax
								weights_to_place[dd,cc] = jeemax
							end
						end
					end #end ltp
				end #end placecell_pos
			end
		end


		###################################################################
		# go over sequence cells
		###################################################################
		# These variables collect new input spikes for the neuron
		forwardInputsEaff[:] = 0.
		forwardInputsErec[:] = 0.
		forwardInputsI[:] = 0.

		# update single cells
		spiked = zeros(Bool,Ncells)

		tprev = dt*(tt-1)
		#over global inhibitory cells first
		for dd = 1:global_inhibitory_cells
			while(t > nextin[dd])
				forwardInputsIPrev[global_inhibitory_excitatory_targets[dd]] += jextin
				forwardInputsIPrev[global_inhibitory_inhibitory_targets[dd]] += jextin * jii/jei0
				if jextin > 0.01
					ns_global_inhibitory_cells[dd] += 1
					times_global_inhibitory_cells[dd, ns_global_inhibitory_cells[dd]] = nextin[dd]
					nextin[dd] += rand(expdist)/rextin
				else
					nextin[dd] += rand(expdist)/default_rextin
				end
			end
		end
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
			for i in 1:size(ext_in_rates)[1]
				if ext_in_rates[i,1] == t
					if cc == 1
						println("ext",ext_in_rates[i,1])
						println("ext2",ext_in_rates[i,2])
						println("ext3",ext_in_rates[i,3])
						println("ext4",ext_in_rates[i,4])
					end
					rextin = ext_in_rates[i,3]
					jextin = ext_in_rates[i,4]
					inhibitory_stimulus_off = true
				end
			end



			#############################################artificially turns on stopping inhibition
			if(inhibitory_stimulus_off)
				if (auxiliary_pattern_recovery_count[auxiliary_pattern_index] < 2)
					if(tt > auxiliary_time[auxiliary_pattern_index] + 120) ### assumes inhibitory stimulus is turned off only after 12ms
						if (cc == 1 && spike_count_over_time_per_pattern[auxiliary_pattern_index+1,tt] > 5) #this number might have to be adapted with the system
							println("spike_count_over_time_per_pattern[auxiliary_pattern_index+1,tt] = ",spike_count_over_time_per_pattern[auxiliary_pattern_index+1,tt])
						end
						if (spike_count_over_time_per_pattern[auxiliary_pattern_index+1,tt] > 10) #this number might have to be adapted with the system
							auxiliary_pattern_recovery_count[auxiliary_pattern_index] += 1
							if(auxiliary_pattern_recovery_count[auxiliary_pattern_index] == 1)
								auxiliary_time[auxiliary_pattern_index] = tt
							end
							if(auxiliary_pattern_recovery_count[auxiliary_pattern_index] == 2)
								auxiliary_pattern_index += 1
								turn_on_time[auxiliary_pattern_index] = t;
								rextin = default_rextin
								jextin = default_jextin
								inhibitory_stimulus_off = false
								println(auxiliary_pattern_index)
								#auxiliary_pattern_index += 1
								println("inh on auxiliary_pattern_index = ",auxiliary_pattern_index)
								println("auxiliary_pattern_index_place = ",auxiliary_pattern_index_place)
								if (auxiliary_pattern_index_place == auxiliary_pattern_index)
									for index = 1:size(cellpathsequence[:,index_cellpathsequence],1) ### should change here if using more than one place cell per position
										for placecell_pos_index = 1:placecells_per_position 
											ee = convert(Int64,(placecell_pos_index-1)*N_placecells/placecells_per_position) + cellpathsequence[index,index_cellpathsequence]
											#cc = cellpathsequence[index,index_cellpathsequence]
											# Next external spike
											rx_place[ee] = rex_place
											nextx_place[ee] = t + rand(expdist)/rex_place
										end

									end
								end
							end
						end
					end
				end
			end
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
            		if cc <= Ne/number_patterns  || Ne+1 <= cc <= Ne + Ni/(number_patterns+1)
               			while(t > nextxstim[cc]) 
                    			nextxstim[cc] += rand(expdist)/stim_rate
                    			if stim_weight > 0
			                        forwardInputsEPrevrec[cc] += stim_weight
			                else
                        			forwardInputsIPrev[cc] -= stim_weight
			                end
		                end
 			end



			# plasticity variables
			if cc <= Ne
				u_vstdp[cc] += dt*(v[cc] - u_vstdp[cc])/tauu;
				v_vstdp[cc] += dt*(v[cc] - v_vstdp[cc])/tauv;
				x_vstdp[cc] -= dt*x_vstdp[cc]/taux;
			end

			# Synaptic rise/decay update
			xeriseaff[cc] += -dt*xeriseaff[cc]/tauerise + forwardInputsEPrevaff[cc]
			xedecayaff[cc] += -dt*xedecayaff[cc]/tauedecay + forwardInputsEPrevaff[cc]
			# Synaptic rise/decay update
			xeriserec[cc] += -dt*xeriserec[cc]/tauerise + forwardInputsEPrevrec[cc]
			xedecayrec[cc] += -dt*xedecayrec[cc]/tauedecay + forwardInputsEPrevrec[cc]
			
			xirise[cc] += -dt*xirise[cc]/tauirise + forwardInputsIPrev[cc]
			xidecay[cc] += -dt*xidecay[cc]/tauidecay + forwardInputsIPrev[cc]


            # Adaptation current for excitatory neurons
			if cc <= Ne
				vth[cc] += dt*(vth0 - vth[cc])/tauth;
				wadapt[cc] += dt*(aw_adapt*(v[cc]-vleake) - wadapt[cc])/tauw_adapt;
			end

			geaff = (xedecayaff[cc] - xeriseaff[cc])/(tauedecay - tauerise);
			geaff_conductance_for_nmda_spikes[cc, 1:dendriticintegrationwindow] = geaff_conductance_for_nmda_spikes[cc, 2:(dendriticintegrationwindow+1)]
			geaff_conductance_for_nmda_spikes[cc,dendriticintegrationwindow+1] = geaff
			gerec = (xedecayrec[cc] - xeriserec[cc])/(tauedecay - tauerise);
			gi = (xidecay[cc] - xirise[cc])/(tauidecay - tauirise);


            # Not in refractory period
			if t > (lastSpike[cc] + taurefrac) 
				# Update excitatory and inhibitory conductance
				geaff = (xedecayaff[cc] - xeriseaff[cc])/(tauedecay - tauerise);

				if (tt > dendriticspiketimes[cc]+dendriticrefractoryperiod+dendriticlatency+30 && tt > dendriticintegrationwindow)
					geaff_slope = geaff_conductance_for_nmda_spikes[cc,dendriticintegrationwindow+1] - geaff_conductance_for_nmda_spikes[cc,1]
					if (geaff_slope >= nmdaspikethreshold)
						dendriticspiketimes[cc] = tt
					end
				end
				if (dendriticspiketimes[cc]+dendriticlatency<=tt<dendriticspiketimes[cc]+30+dendriticlatency)#30 is the length of the dendritic spike
					dendritic_time_index = convert(Int64,tt-dendriticspiketimes[cc]-dendriticlatency) + 1
					geaff = 4*(10^3)*dendriticInputcurve[dendritic_time_index]/(erev-v[cc])

				end 
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

					if ns[cc] <= Nspikes
						times[cc,ns[cc]] = t;
					end

					if cc<= Ne
						x_vstdp[cc] += 1./taux;
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

					# Loop over synaptic projections - sequence cells # unnecessary
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
					
					# Loop over synaptic projections - place cells
					for dd = 1:N_placecells
						if cc <= Ne # Excitatory synapse
							forwardInputsEplace[dd] += weights_to_place[cc,dd] 
						end
					end
				end #end if(spiked)
			end #end if(not refractory)

			##########################################################################################
			# plasticity sequence to place cells
			##########################################################################################
			#vstdp, ltd component
			if spiked[cc] && (cc <= Ne)
				for dd = 1:N_placecells #depress weights from cc to cj
					if weights_to_place[cc,dd] == 0.
						continue
					end
					if u_vstdp_place[dd] > thetaltd
						weights_to_place[cc,dd] -= altd*(u_vstdp_place[dd]-thetaltd)
						if weights_to_place[cc,dd] < jeemin
							weights_to_place[cc,dd] = jeemin
						end
					end
				end
			end #end ltd

		end #end loop over cells
		

		forwardInputsEPrevaff = copy(forwardInputsEaff)
		forwardInputsEPrevrec = copy(forwardInputsErec)
		forwardInputsIPrev = copy(forwardInputsI)
		forwardInputsEPrevplace = copy(forwardInputsEplace)

	end #end loop over time
	@printf("\r")


	if maximum(ns) > Nspikes
		println("Warning! Maximum number of spikes observed was ",maximum(ns)," whereas memory was only allocated for ",Nspikes," spikes.")
	end
	times = times[:,1:maximum(ns)]
	times_place = times_place[:,1:maximum(ns_place)]

	return times,ns,time_xaxis, spike_count_over_time_per_pattern,turn_off_times, turn_on_time,ex_spike_count_per_step_per_pattern,in_spike_count_per_step_per_pattern, times_place,ns_place,time_xaxis,cellpathsequence_times, weights_to_place, times_global_inhibitory_cells, ns_global_inhibitory_cells, weights
end


