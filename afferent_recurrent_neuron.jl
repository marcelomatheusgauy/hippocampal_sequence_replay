# This file is made by Marcelo Gauy
# this function receives two neuron indices cc and dd, checks whether they are excitatory and their corresponding pattern indices
# if the pattern indices are consecutive then there could be an afferent synapse between the 2 neurons so the function returns true.
# all other cases it returns false

function is_aff(cc,dd,Ne,number_patterns)
	cc_pattern_number = -50;
	dd_pattern_number = -100;
	if (cc <= Ne*number_patterns)
		cc_pattern_number = floor(cc/Ne) + 1;
	end

	if (dd <= Ne*number_patterns)
		dd_pattern_number = floor(dd/Ne) + 1;
	end

	if (dd_pattern_number == cc_pattern_number + 1 || dd_pattern_number == cc_pattern_number - 1)
		return true
	else
		return false
	end
end
