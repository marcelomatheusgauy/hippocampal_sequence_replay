# This file is made by Marcelo Gauy

# given the spiking time of every cell, computes the number of spikes present for a cell in a given time frame;
# the returned structure is indexed by place cell position!

function spike_count(times,framelength, totalsimlength,gridsize, Ncells, placecells_per_position)

	secondcoordlength = convert(Int64,floor(totalsimlength/framelength))
	cell_spike_count_over_frame = zeros(gridsize,gridsize,secondcoordlength)
	spike_count_per_cell = zeros(Int64,Ncells,1)

	for index = 1:secondcoordlength
		for cc = 1:Ncells
			place_pos_index = cc - convert(Int64,floor((cc-1)/(gridsize*gridsize)))*gridsize*gridsize
			x_coord = convert(Int64,floor((place_pos_index-1)/gridsize)+1)#div(cc,gridsize)
			y_coord = place_pos_index%gridsize
			if (y_coord == 0)
				y_coord = gridsize
			end
			spike = spike_count_per_cell[cc]+1
			max_spikes = size(times,2)

			this_frame = true

			while(this_frame)
				if (spike <= max_spikes)
					if((index-1)*framelength < times[cc,spike] <= index*framelength)
						cell_spike_count_over_frame[x_coord,y_coord,index] += 1
						spike += 1
					else
						this_frame = false
					end
				else
					this_frame = false
				end
				spike_count_per_cell[cc] = spike - 1
			end
		end
	end
	
	return cell_spike_count_over_frame

end
