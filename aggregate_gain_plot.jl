# This file is made by Marcelo Gauy and is based on
# simulation code from Litwin-Kumar 2014.

# In this simulation we model two external rates!
#   First corresponds to the previous external rate
#   Second corresponds to the excitatory population input whose rate we vary
# The blue curve corresponds to excitatory gain and the red to inhibitory gain.
# The black curve is the identity for reference since we want the excitatory gain
# to be below the identity initially and then cross it so it can self-sustain.

# the aim is to build the phase plane of the excitatory and inhibitory population 
# rates to analyse how the bistable unit behaves. This is first done without any external input
# This is adapted from the gainfunction.jl

using PyPlot
using PyCall

@pyimport numpy

include("simulation7.jl")
include("parameters.jl")


function weightpars_custom() #parameters needed to generate weight matrix
	Ne = 1 # Number of excitatory neurons
	Ni = 1 # Number of inhibitory neurons
	jee0 = 10.0 # initial ee strength
	jei0 = 10.0 # inhibition to excitation strength
	jie = 4.0 # excitation to inhibition strength
	jii = 3.0 # ii strength 
	#jba = 0.9 # A to B strength
	p = 0.2
	rex = 1.0 # external input rate to e (khz) 
	rix = 1.0 # external input rate to i (khz)
	jex = 2.0  # external to e strength
	jix = -4.0 # external to i strength
	rextin = 1.0 #stop inhibition
	jextin = 2.5
	rate_neighbour_pattern = 1.8 #can be used to simulate the effect of the incoming weight of the previous pattern, or that of the pattern to come
	weight_neighbour_pattern = 1.6 
	return Ne,Ni,jee0,jei0,jie,jii,p,rex,rix,jex,jix,rextin,jextin, rate_neighbour_pattern, weight_neighbour_pattern
end

Ne,Ni,_,_,_,_,p,_,_,_,_ = weightpars_custom()

function stimpars_custom()
    return [1.0 24.0 0.25 0.0]
end

_,_,stim_rate = stimpars_custom()

function simpars_custom()
    #simulation
	dt = .1 #integration timestep
	T = 2000 #simulation time
	Nskip = 1000 #how often (in number of timesteps) to save w_in
	vpeak = 20 #cutoff for voltage.  when crossed, record a spike and reset
	Nspikes = 10000 #maximum number of spikes to record per neuron

    return dt, T, Nskip, vpeak, Nspikes
end

_,T,_,_,_ = simpars_custom()

# The number of excitatory neurons we simulate
nex = 200
nin = 100

tic()
datestring = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")

ex_rates = linspace(0.0,0.25,401)
in_rates = linspace(0.0,0.25,401)

Y1, Y2 = numpy.meshgrid(1000*ex_rates,1000*in_rates)

rates_ex = zeros(length(ex_rates),length(in_rates))
rates_in = zeros(length(ex_rates),length(in_rates))
rates_ex_diff = zeros(length(ex_rates),length(in_rates))
rates_in_diff = zeros(length(ex_rates),length(in_rates))

diagonal = [ex_rates, ex_rates]

i = 1
##look for implicit solution - that is the point of inhibitory rate i which solves i = g(e,i)
implicit_inhibitory_rate_stable = zeros(length(ex_rates))
implicit_excitatory_rate_stable = zeros(length(ex_rates))

for ex_rate in ex_rates
	j = 1
	println(i)
	for in_rate in in_rates
		function internalpars_custom()
			internal_input = true
			r_int = nex*p*ex_rate
			ri_int = nin*p*in_rate
			return internal_input,r_int,ri_int
		end
    		times,ns = sim(weightpars_custom,internalpars_custom,membranepars,synapsepars,stimpars_custom,simpars_custom)
		#println("Excitatory firing rate: ",1000*ns[1]/T," Hz")
		#println("Inhibitory firing rate: ",1000*ns[2]/T," Hz")

		rates_ex[i,j] = 1000*ns[1]/T
		rates_in[i,j] = 1000*ns[2]/T

		if rates_in[i,j] - 1000*in_rates[j] > -4.0 && rates_in[i,j] - 1000*in_rates[j] < 4.0 && rates_in[i,j] > 8.0
			implicit_inhibitory_rate_stable[i] = in_rates[j]
		end

		if rates_ex[i,j] - 1000*ex_rates[i] > -4.0 && rates_ex[i,j] - 1000*ex_rates[i] < 4.0 && rates_ex[i,j] > 8.0
			implicit_excitatory_rate_stable[j] = ex_rates[i]
			#println("implicit_excitatory_rate_stable[i] = ",implicit_excitatory_rate_stable[i])
		end

		j += 1
	end
	i += 1
end

##Find excitatory gain function given implicit inhibitory rate
excitatory_gain_function = zeros(length(ex_rates))

i = 1
for ex_rate in ex_rates
	function internalpars_custom()
		internal_input = true
		r_int = nex*p*ex_rate
		ri_int = nin*p*implicit_inhibitory_rate_stable[i]
		return internal_input,r_int,ri_int
	end
	times,ns = sim(weightpars_custom,internalpars_custom,membranepars,synapsepars,stimpars_custom,simpars_custom)
	println("Excitatory firing rate: ",1000*ns[1]/T," Hz")
	println("Inhibitory firing rate: ",1000*ns[2]/T," Hz")

	excitatory_gain_function[i] = 1000*ns[1]/T

	i += 1
end


##Find inhibitory gain function given implicit inhibitory rate
inhibitory_gain_function = zeros(length(ex_rates))

i = 1
for in_rate in in_rates
	function internalpars_custom()
		internal_input = true
		r_int = nex*p*implicit_excitatory_rate_stable[i]
		ri_int = nin*p*in_rate
		return internal_input,r_int,ri_int
	end
	times,ns = sim(weightpars_custom,internalpars_custom,membranepars,synapsepars,stimpars_custom,simpars_custom)
	println("Excitatory firing rate: ",1000*ns[1]/T," Hz")
	println("Inhibitory firing rate: ",1000*ns[2]/T," Hz")


	inhibitory_gain_function[i] = 1000*ns[2]/T

	i += 1
end

figure(figsize=(4,4))
plot(1000*ex_rates,excitatory_gain_function, color="red", linewidth=2.0, linestyle="-")
plot(1000*ex_rates,inhibitory_gain_function, color="teal", linewidth=2.0, linestyle="-")
#plot(1000*ex_rates,1000*implicit_excitatory_rate_stable, color="blue", linewidth=2.0, linestyle="-")
plot(1000*ex_rates, 1000*ex_rates, color="black", linewidth=2.0, linestyle="-")
xlim(0,70)
ylim(0,70)
ylabel("Rate (Hz)", fontsize = 17)
xlabel("Rate (Hz)", fontsize = 17)
xticks(fontsize = 15)
yticks(fontsize = 15)
tight_layout()
savefig(string("figures/",datestring,"aggregate_gain_function.png"),dpi=150)

toc()


