# This file is made by Hafsteinn Einarsson and is based on
# simulation code from Litwin-Kumar 2014.

# This file contains default parameter values, please do not change the constants in this file.
# In order to change parameters create a new function and pass it to simulation.jl. For an example
# of this see gainfunction.jl or balanced_activity_raster.jl.

# Please note that Litwin-Kumar uses the naming
# convention that jie means weight from excitatory
# to inhibitory...
function weightpars() #parameters needed to generate weight matrix
	Ne = 4000 # Number of excitatory neurons
	Ni = 1000 # Number of inhibitory neurons
	jee0 = 2.86 #initial ee strength
	jei0 = 48.7 #initial ei strength
	jie = 1.27 #ie strength (not plastic)
	jii = 16.2 #ii strength (not plastic)
	p = 0.2

	rex = 4.5 #external input rate to e (khz)
	rix = 2.25 #external input rate to i (khz)

	jex = 1.78 #external to e strength
	jix = 1.27 #external to i strength
	return Ne,Ni,jee0,jei0,jie,jii,p,rex,rix,jex,jix
end

# This simulates input from an excitatory population
# at rate r_int. The weights used are jee0 and jie.
# If using an excitatory population this is not necessary
# since we also have external input. This is just used
# at the moment to compute the gain function.
function internalpars()
	internal_input = false
	r_int = 0.0
	ri_int = 0.0
	return internal_input,r_int,ri_int
end

function membranepars()
	#membrane dynamics
	taue = 20 #e membrane time constant
	taui = 20 #i membrane time constant
	vleake = -70 #e resting potential
	vleaki = -62 #i resting potential
	deltathe = 2 #eif slope parameter
	C = 300 #capacitance
	erev = 0 #e synapse reversal potential
	irev = -75 #i synapse reversal potntial
	vth0 = -52 #initial spike voltage threshold
	ath = 10 #increase in threshold post spike
	tauth = 30 #threshold decay timescale
	vre = -60 #reset potential
	taurefrac = 1 #absolute refractory period originally 1
	aw_adapt = 4 #adaptation parameter a
	bw_adapt = .805 #adaptation parameter b #originally .805
	tauw_adapt = 150 #adaptation timescale
	return taue,taui,vleake,vleaki,deltathe,C,erev,irev,vth0,ath,tauth,vre,taurefrac,aw_adapt,bw_adapt,tauw_adapt
end

function synapsepars()
	tauerise = 1 #e synapse rise time it was 1 in Litwin-Kumar
	tauedecay = 6 #e synapse decay time it was 6 in Litwin-Kumar
	tauirise = .5 #i synapse rise time
	tauidecay = 2 #i synapse decay time #it was originally 2 in Litwin-Kumar
	nmdaspikethreshold = 6.0 #originally 6.0 -> #increase this to avoid dendritic spikes. Used in Figure 3
	####dendritic spikes
	dendriticA = 55 #nA
	dendriticB = 64 #nA
	dendriticC = 9 #nA
	taudendriticA = 2 #ms/10
	taudendriticB = 3 #ms/10
	taudendriticC = 7 #ms/10
	dendriticlatency = 27 #ms/10 it is 27 in Jahnke
	dendriticintegrationwindow = 20 #ms/10
	dendriticrefractoryperiod = 50 #ms/10
	return tauerise,tauedecay,tauirise,tauidecay,nmdaspikethreshold,dendriticA,dendriticB,dendriticC,taudendriticA,taudendriticB,taudendriticC,dendriticlatency, dendriticintegrationwindow, dendriticrefractoryperiod
end

# Columns are: (1) Rate, (2) start time, (3) end time
# Current implementation does not allow overlapping stimuli
function stimpars()
	return [] # No stimulus by default
end

function simpars()
    #simulation
	dt = .1 #integration timestep
	T = 2000 #simulation time
	Nskip = 1000 #how often (in number of timesteps) to save w_in
	vpeak = 20 #cutoff for voltage.  when crossed, record a spike and reset
	Nspikes = 10000 #maximum number of spikes to record per neuron

    return dt, T, Nskip, vpeak, Nspikes
end

function voltagestdppars()
	altd = .0008 #ltd strength
	altp = .0014 #ltp strength
	thetaltd = -70 #ltd voltage threshold
	thetaltp = -49 #ltp voltage threshold
	tauu = 10 #timescale for u variable
	tauv = 7 #timescale for v variable
	taux = 15 #timescale for x variable
	return altd,altp,thetaltd,thetaltp,tauu,tauv,taux
end

function inhstdppars()
	tauy = 20 #width of istdp curve
	eta = 1 #istdp learning rate
	r0 = .003 #target rate (khz)
	return tauy,eta,r0
end
