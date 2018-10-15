# hippocampal_sequence_replay
A hippocampal model for behavioral time acquisition and fast bidirectional replay of spatio-temporal memory sequences

This github project contains the codes used in the simulations of the paper: A hippocampal model for behavioral time acquisition and fast bidirectional replay of spatio-temporal memory sequences. We make the codes public, so the results obtained can be reproduced. In the following, I explain the file structure:

sequencereplayencodingmode_fast_replay_animation_linear_track_loops.jl: To simulate the entire model just run this file(julia); this simulates the model (encoding including place and sequence cells) with the parameters in parameters.jl (and those provided within the file itself); after the encoding mode is finished, it simulates a single replay and reverse replay of place and sequence cells; the data is saved (in the data folder) and preliminary plots for encoding, replay and reverse replay are produced (in figures folder). We note that the data file produced is too large and we do not include a preliminary one in the data folder.

parameters.jl: contains the basic neuronal parameters used. These are changed externally on convenience in the relevant files.

simulationsequencereplayencodingmode.jl: contains the function that actually implements the encoding mode of the model. Is called by sequencereplayencodingmode_fast_replay_animation_linear_track_loops.jl; Note that the random seed is fixed to ensure reproducibility.

simulationfastreplaymode.jl: contains the function that actually implements the forward replay mode of the model. Is called by sequencereplayencodingmode_fast_replay_animation_linear_track_loops.jl; Note that the random seed is fixed to ensure reproducibility.

simulationreversereplay.jl: contains the function that actually implements the reverse replay of the model. Is called by sequencereplayencodingmode_fast_replay_animation_linear_track_loops.jl; Note that the random seed is fixed to ensure reproducibility.

aggregate_spike_statistics_over_time.jl, afferent_recurrent_neuron.jl and mousepathtocellpathsequence.jl are specific utility functions called when necessary.

loadscript_linear_loops.jl contains an example code which produces plots that look like the ones in the paper (Figure 1). Note that we do not provide the data file for the figures attached because it is too large; To get the data file, one must run sequencereplayencodingmode_fast_replay_animation_linear_track_loops.jl.

In case of doubts please write to: marcelo.matheus@inf.ethz.ch
