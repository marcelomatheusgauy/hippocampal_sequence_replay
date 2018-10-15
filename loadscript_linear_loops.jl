using HDF5

using IJulia
using PyCall
@pyimport matplotlib.animation as anim
using PyPlot


using JLD

d = load("data/video_data_4.jld");

times = d["times"];
ns = d["ns"];
times_x_axis = d["times_x_axis"];
times_place = d["times_place"];
ns_place = d["ns_place"];
cellpathsequence_times = d["cellpathsequence_times"];
weights_to_place = d["weights_to_place"];
mousepathpositions = d["mousepathpositions"];
cellpathsequence = d["cellpathsequence"];
times_fast = d["times_fast"];
ns_fast = d["ns_fast"];
times_place_fast = d["times_place_fast"];
ns_place_fast = d["ns_place_fast"];
times_reverse = d["times_reverse"];
ns_reverse = d["ns_reverse"];
times_place_reverse = d["times_place_reverse"];
ns_place_reverse = d["ns_place_reverse"];
turn_on_time = d["turn_on_time"];
mousemovingtimes = d["mousemovingtimes"];
mousestoppingtimes = d["mousestoppingtimes"];
linear_track_size = d["linear_track_size"];


N_placecells = size(ns_place,1)
Nsteps = size(times_x_axis,1)
T = convert(Int64,Nsteps*0.1)
placecells_per_position = 100

println(linear_track_size)
println(N_placecells)
println(Nsteps)
println(T)


#####spike count
spike_count_place_per_time_per_pos = zeros(linear_track_size, T)
spike_count_place_per_cell = zeros(Int,N_placecells)
for time_index = 1:T
	for pos_index = 1:linear_track_size
		for place_pos_index = 1:placecells_per_position
			cc = convert(Int64,(place_pos_index-1)*N_placecells/placecells_per_position) + pos_index
			next_spike = spike_count_place_per_cell[cc] + 1
			if next_spike <= ns_place[cc]
				if 0.0 < times_place[cc,next_spike] <= time_index
					spike_count_place_per_time_per_pos[pos_index,time_index] += 1.0
					spike_count_place_per_cell[cc] += 1
				end
			end
		end
	end
end
###################################
#####smoothing out the spike count
smoothing_time = 5 #divisor of T
place_rate_over_time_smoothed = zeros(linear_track_size, T)
for time_index = 1:T
	if (time_index <= smoothing_time)
		aux_spike_count = sum(spike_count_place_per_time_per_pos[:, 1:(time_index+smoothing_time)],2) # the 10 is 1/dt of the simulation
		place_rate_over_time_smoothed[:,time_index] = 1000*aux_spike_count/(placecells_per_position*(time_index+smoothing_time+1))
	elseif time_index < T - smoothing_time
		aux_spike_count = sum(spike_count_place_per_time_per_pos[:, (time_index-smoothing_time):(time_index+smoothing_time)],2)
		place_rate_over_time_smoothed[:,time_index] = 1000*aux_spike_count/(placecells_per_position*(2*smoothing_time+1))
	else
		aux_spike_count = sum(spike_count_place_per_time_per_pos[:, (time_index-smoothing_time):T],2) # the 10 is 1/dt of the simulation
		place_rate_over_time_smoothed[:,time_index] = 1000*aux_spike_count/(placecells_per_position*(T-time_index+smoothing_time+1))
	end
end
###################################

track_size = 12
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
ax[:set_xlabel]("Time(s)", fontsize=35)
ax[:set_ylabel]("Track Position", fontsize=35)
axsec = ax[:twinx]()
axsec[:set_ylabel]("Rate(Hz)", fontsize=35)
axsec[:spines]["top"][:set_color]("none")
axsec[:spines]["bottom"][:set_color]("none")
axsec[:spines]["left"][:set_color]("none")
axsec[:spines]["right"][:set_color]("none")
axsec[:tick_params](labelcolor="w", top="off", bottom="off", left="off", right="off")
for pos_index = 1:track_size
	##excitatory raster plots
	rowcount = 1
	if (pos_index < 13)
		fig[:add_subplot](track_size, 1, pos_index, sharex = ax)
		ax1 = gca()
		if (pos_index< track_size)
			setp(ax1[:get_xticklabels](),visible=false) # Disable x tick labels
		else
			ax1[:tick_params](labelsize=20)
			ax1[:set_xticklabels](["0","0","2","4","6","8","10","12","14"])
		end
		#for index_cc = 1:placecells_per_position
		#	cc = convert(Int64,(index_cc-1)*N_placecells/placecells_per_position) + pos_index
		#	vals = times_place[cc,1:ns_place[cc]]
		#	y = rowcount*ones(length(vals))
		#	scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
			#ax[:scatter](vals,y,s=.3,c="k",marker="o",linewidths=0)
		#	rowcount+=1
		#end
	elseif (pos_index > 13)
		fig[:add_subplot](track_size, 1, pos_index-1, sharex = ax)
		ax1 = gca()
		if (pos_index< track_size)
			setp(ax1[:get_xticklabels](),visible=false) # Disable x tick labels
		else
			ax1[:tick_params](labelsize=20)
			ax1[:set_xticklabels](["0","0","2","4","6","8","10","12","14"])
		end
		#for index_cc = 1:placecells_per_position
		#	cc = convert(Int64,(index_cc-1)*N_placecells/placecells_per_position) + pos_index
		#	vals = times_place[cc,1:ns_place[cc]]
		#	y = rowcount*ones(length(vals))
		#	scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
			#ax[:scatter](vals,y,s=.3,c="k",marker="o",linewidths=0)
		#	rowcount+=1
		#end
	end
	if pos_index != 13
		xlim(0,T)
		setp(ax1[:get_yticklabels](),visible=false) # Disable x tick labels
		ylim(0,placecells_per_position)
		if pos_index == 1
			ylabel("\$A\$", fontsize=25,fontdict=Dict("color"=>"teal"))
		elseif pos_index == 2
			ylabel("\$B\$", fontsize=25,fontdict=Dict("color"=>"brown"))
		elseif pos_index == 3
			ylabel("\$C\$", fontsize=25,fontdict=Dict("color"=>"red"))
		elseif pos_index == 4
			ylabel("\$D\$", fontsize=25,fontdict=Dict("color"=>"blue"))
		elseif pos_index == 5
			ylabel("\$E\$", fontsize=25)
		elseif pos_index == 6
			ylabel("\$F\$", fontsize=25)
		elseif pos_index == 7
			ylabel("\$G\$", fontsize=25)
		elseif pos_index == 8
			ylabel("\$H\$", fontsize=25)
		elseif pos_index == 9
			ylabel("\$I\$", fontsize=25)
		elseif pos_index == 10
			ylabel("\$J\$", fontsize=25)
		elseif pos_index == 11
			ylabel("\$K\$", fontsize=25)
		elseif pos_index == 12
			ylabel("\$L\$", fontsize=25)
		elseif pos_index == 13
			ylabel("\$M\$", fontsize=25)
		elseif pos_index == 14
			ylabel("\$N\$", fontsize=25)
		elseif pos_index == 15
			ylabel("\$14\$", fontsize=25)
		elseif pos_index == 16
			ylabel("\$15\$", fontsize=25)
		end
		if pos_index == 1
			ax2 = ax1[:twinx]()
			ax2[:plot](place_rate_over_time_smoothed[pos_index,:],c="teal")
			#ax2[:set_ylabel]("Rate", color = "b")
			#setp(ax2[:get_yticklabels](),visible=false) # Disable x tick labels
			ax2[:tick_params](colors="teal", labelsize = 15)
		elseif pos_index == 2
			ax2 = ax1[:twinx]()
			ax2[:plot](place_rate_over_time_smoothed[pos_index,:],c="brown")
			#ax2[:set_ylabel]("Rate", color = "b")
			#setp(ax2[:get_yticklabels](),visible=false) # Disable x tick labels
			ax2[:tick_params](colors="brown", labelsize = 15)
		elseif pos_index == 3
			ax2 = ax1[:twinx]()
			ax2[:plot](place_rate_over_time_smoothed[pos_index,:],c="r")
			#ax2[:set_ylabel]("Rate", color = "b")
			#setp(ax2[:get_yticklabels](),visible=false) # Disable x tick labels
			ax2[:tick_params](colors="r", labelsize = 15)
		elseif pos_index == 4
			ax2 = ax1[:twinx]()
			ax2[:plot](place_rate_over_time_smoothed[pos_index,:],c="blue")
			#ax2[:set_ylabel]("Rate", color = "b")
			#setp(ax2[:get_yticklabels](),visible=false) # Disable x tick labels
			ax2[:tick_params](colors="blue", labelsize = 15)
		else
			ax2 = ax1[:twinx]()
			ax2[:plot](place_rate_over_time_smoothed[pos_index,:], c="k")
			#ax2[:set_ylabel]("Rate", color = "b")
			#setp(ax2[:get_yticklabels](),visible=false) # Disable x tick labels
			ax2[:tick_params](labelsize = 15)
		end
		ylim(0,80)
	end
end
tight_layout()
subplots_adjust(hspace=0.0)


savefig(string("figures/",datestring,"_raster_for_place_linear_track.png"),dpi=300)
#


T_fast = 100

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
ax[:set_xlabel]("Time(ms)",fontsize=35)
ax[:set_ylabel]("Track Position",fontsize=35)
#axsec = ax[:twinx]()
#axsec[:set_ylabel]("Rate(Hz)", color = "b",fontsize=35)
#axsec[:spines]["top"][:set_color]("none")
#axsec[:spines]["bottom"][:set_color]("none")
#axsec[:spines]["left"][:set_color]("none")
#axsec[:spines]["right"][:set_color]("none")
#axsec[:tick_params](labelcolor="w", top="off", bottom="off", left="off", right="off")
for pos_index = 1:track_size
	##excitatory raster plots
	rowcount = 1
	if (pos_index < 13)
		fig[:add_subplot](track_size, 1, pos_index, sharex = ax)
		ax1 = gca()
		if (pos_index< track_size)
			setp(ax1[:get_xticklabels](),visible=false) # Disable x tick labels
		else
			ax1[:tick_params](labelsize=20)
		end
		for index_cc = 1:placecells_per_position
			cc = convert(Int64,(index_cc-1)*N_placecells/placecells_per_position) + pos_index
			vals = times_place_fast[cc,1:ns_place_fast[cc]]
			y = rowcount*ones(length(vals))
			if pos_index == 1
				scatter(vals,y,s=2,c="teal",marker="o",linewidths=0)
			elseif pos_index == 2
				scatter(vals,y,s=2,c="brown",marker="o",linewidths=0)
			elseif pos_index == 3
				scatter(vals,y,s=2,c="r",marker="o",linewidths=0)
			elseif pos_index == 4
				scatter(vals,y,s=2,c="b",marker="o",linewidths=0)
			else
				scatter(vals,y,s=2,c="k",marker="o",linewidths=0)
			end
			#ax[:scatter](vals,y,s=.3,c="k",marker="o",linewidths=0)
			rowcount+=1
		end
	elseif (pos_index > 13)
		fig[:add_subplot](track_size, 1, pos_index-1, sharex = ax)
		ax1 = gca()
		if (pos_index< track_size)
			setp(ax1[:get_xticklabels](),visible=false) # Disable x tick labels
		else
			ax1[:tick_params](labelsize=20)
		end
		for index_cc = 1:placecells_per_position
			cc = convert(Int64,(index_cc-1)*N_placecells/placecells_per_position) + pos_index
			vals = times_place_fast[cc,1:ns_place_fast[cc]]
			y = rowcount*ones(length(vals))
			scatter(vals,y,s=2,c="k",marker="o",linewidths=0)
			#ax[:scatter](vals,y,s=.3,c="k",marker="o",linewidths=0)
			rowcount+=1
		end
	end
	if pos_index != 13
		xlim(0,T_fast)
		setp(ax1[:get_yticklabels](),visible=false) # Disable x tick labels
		ylim(0,placecells_per_position)
		if pos_index == 1
			ylabel("\$A\$", fontsize=25,fontdict=Dict("color"=>"teal"))
		elseif pos_index == 2
			ylabel("\$B\$", fontsize=25,fontdict=Dict("color"=>"brown"))
		elseif pos_index == 3
			ylabel("\$C\$", fontsize=25,fontdict=Dict("color"=>"red"))
		elseif pos_index == 4
			ylabel("\$D\$", fontsize=25,fontdict=Dict("color"=>"blue"))
		elseif pos_index == 5
			ylabel("\$E\$", fontsize=25)
		elseif pos_index == 6
			ylabel("\$F\$", fontsize=25)
		elseif pos_index == 7
			ylabel("\$G\$", fontsize=25)
		elseif pos_index == 8
			ylabel("\$H\$", fontsize=25)
		elseif pos_index == 9
			ylabel("\$I\$", fontsize=25)
		elseif pos_index == 10
			ylabel("\$J\$", fontsize=25)
		elseif pos_index == 11
			ylabel("\$K\$", fontsize=25)
		elseif pos_index == 12
			ylabel("\$L\$", fontsize=25)
		elseif pos_index == 13
			ylabel("\$M\$", fontsize=25)
		elseif pos_index == 14
			ylabel("\$N\$", fontsize=25)
		elseif pos_index == 15
			ylabel("\$14\$", fontsize=25)
		elseif pos_index == 16
			ylabel("\$15\$", fontsize=25)
		end
		#ax2 = ax1[:twinx]()
		#ax2[:plot](place_rate_over_time_smoothed_fast[pos_index,:],c="b")
		#ax2[:set_ylabel]("Rate", color = "b")
		#setp(ax2[:get_yticklabels](),visible=false) # Disable x tick labels
		#ax2[:tick_params](colors="b", labelsize = 15)
	end
end

tight_layout()
subplots_adjust(hspace=0.0)
savefig(string("figures/",datestring,"_raster_for_place_linear_track_replay.png"),dpi=300)
#


T_reverse = 150

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
ax[:set_xlabel]("Time(ms)",fontsize=35)
ax[:set_ylabel]("Track Position",fontsize=35)
#axsec = ax[:twinx]()
#axsec[:set_ylabel]("Rate(Hz)", color = "b",fontsize=35)
#axsec[:spines]["top"][:set_color]("none")
#axsec[:spines]["bottom"][:set_color]("none")
#axsec[:spines]["left"][:set_color]("none")
#axsec[:spines]["right"][:set_color]("none")
#axsec[:tick_params](labelcolor="w", top="off", bottom="off", left="off", right="off")
for pos_index = 1:track_size
	##excitatory raster plots
	rowcount = 1
	if (pos_index < 13)
		fig[:add_subplot](track_size, 1, pos_index, sharex = ax)
		ax1 = gca()
		if (pos_index< track_size)
			setp(ax1[:get_xticklabels](),visible=false) # Disable x tick labels
		else
			ax1[:tick_params](labelsize=20)
			ax1[:set_xticklabels](["0","20","40","60","80","100"])
		end
		for index_cc = 1:placecells_per_position
			cc = convert(Int64,(index_cc-1)*N_placecells/placecells_per_position) + pos_index
			vals = times_place_reverse[cc,1:ns_place_reverse[cc]]
			y = rowcount*ones(length(vals))
			if pos_index == 1
				scatter(vals,y,s=2,c="teal",marker="o",linewidths=0)
			elseif pos_index == 2
				scatter(vals,y,s=2,c="brown",marker="o",linewidths=0)
			elseif pos_index == 3
				scatter(vals,y,s=2,c="r",marker="o",linewidths=0)
			elseif pos_index == 4
				scatter(vals,y,s=2,c="b",marker="o",linewidths=0)
			else
				scatter(vals,y,s=2,c="k",marker="o",linewidths=0)
			end
			#ax[:scatter](vals,y,s=.3,c="k",marker="o",linewidths=0)
			rowcount+=1
		end
	elseif (pos_index > 13)
		fig[:add_subplot](track_size, 1, pos_index-1, sharex = ax)
		ax1 = gca()
		if (pos_index< track_size)
			setp(ax1[:get_xticklabels](),visible=false) # Disable x tick labels
		else
			ax1[:tick_params](labelsize=20)
		end
		for index_cc = 1:placecells_per_position
			cc = convert(Int64,(index_cc-1)*N_placecells/placecells_per_position) + pos_index
			vals = times_place_reverse[cc,1:ns_place_reverse[cc]]
			y = rowcount*ones(length(vals))
			scatter(vals,y,s=2,c="k",marker="o",linewidths=0)
			#ax[:scatter](vals,y,s=.3,c="k",marker="o",linewidths=0)
			rowcount+=1
		end
	end
	if pos_index != 13
		xlim(50,T_reverse)
		setp(ax1[:get_yticklabels](),visible=false) # Disable y tick labels
		ylim(0,placecells_per_position)
		if pos_index == 1
			ylabel("\$A\$", fontsize=25,fontdict=Dict("color"=>"teal"))
		elseif pos_index == 2
			ylabel("\$B\$", fontsize=25,fontdict=Dict("color"=>"brown"))
		elseif pos_index == 3
			ylabel("\$C\$", fontsize=25,fontdict=Dict("color"=>"red"))
		elseif pos_index == 4
			ylabel("\$D\$", fontsize=25,fontdict=Dict("color"=>"blue"))
		elseif pos_index == 5
			ylabel("\$E\$", fontsize=25)
		elseif pos_index == 6
			ylabel("\$F\$", fontsize=25)
		elseif pos_index == 7
			ylabel("\$G\$", fontsize=25)
		elseif pos_index == 8
			ylabel("\$H\$", fontsize=25)
		elseif pos_index == 9
			ylabel("\$I\$", fontsize=25)
		elseif pos_index == 10
			ylabel("\$J\$", fontsize=25)
		elseif pos_index == 11
			ylabel("\$K\$", fontsize=25)
		elseif pos_index == 12
			ylabel("\$L\$", fontsize=25)
		elseif pos_index == 13
			ylabel("\$M\$", fontsize=25)
		elseif pos_index == 14
			ylabel("\$N\$", fontsize=25)
		elseif pos_index == 15
			ylabel("\$14\$", fontsize=25)
		elseif pos_index == 16
			ylabel("\$15\$", fontsize=25)
		end
	end
end

tight_layout()
subplots_adjust(hspace=0.0)
savefig(string("figures/",datestring,"_raster_for_place_linear_track_reverse_replay.png"),dpi=300)
#


