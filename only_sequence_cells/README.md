# hippocampal_sequence_replay

This subfolder contains the codes without simulating any place cell activity. Usually these are better to make tests, as they are simpler and faster. They are included here for completeness. Below is a brief description on how to run the files:

turn_ext_inhibition_off_trigger_detector_many_patterns.jl: To produce data corresponding to the encoding mode of the model run this file(julia); this simulates the model (encoding) with the parameters in parameters.jl (and those provided within the file itself), saves the data (in data folder) and produces a preliminary plot (in figures folder).

synfiremode.jl: to produce data corresponding to forward replays run this file; this simulates the replay mode of the model, saves the data and produces a preliminary plot.

reversalreplay.jl: to produce data corresponding to reverse replays run this file; this simulates the reverse replay mode of the model, saves the data and produces a preliminary plot.

parameters.jl: contains the basic neuronal parameters used. These are changed externally on convenience in the relevant files; To remove dendritic spikes, just increase the parameter nmdaspikethreshold from 6.0 to any high value (say 600.0)

simulation5.jl: contains the function that actually implements the encoding mode of the model. Is called by turn_ext_inhibition_off_trigger_detector_many_patterns.jl;

simulation6.jl: contains the function that actually implements the forward replay mode of the model. Is called by synfiremode.jl;

simulation8.jl: contains the function that actually implements the reverse replay of the model. Is called by reversalreplay.jl;

aggregate_spike_statistics_over_time.jl and afferent_recurrent_neuron.jl are specific utility functions called when necessary.

The load_scripts contain example codes the produce plots which look like some of the ones used in the paper.

In case of doubts please write to: marcelo.matheus@inf.ethz.ch
