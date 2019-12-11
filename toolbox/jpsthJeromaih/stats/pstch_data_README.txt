pstch_vector 			= n x m matrix, where n is the number of pairs in that session and m is 1200 and 
						NaN-padded. Each column is a pstch value starting at -50ms and ending at 
						90% of the longest SRT.

pstch_area_all			= the area under the pstch from the time of target onset to the 90% longest SRT, 
						normalized for the length of time it spans. 

no_tdt_pstch_area 		= if this pair is one in which neither cell selects the target, then this is the 
						same value as pstch_area_all for that pair
						
area_before_tdt_xor		= if only one cell selects the target, this is the pstch area between target onset and 
						TDT, normalized by TDT
						
area_before_tdt_both 	= if both cells select the target, this is the pstch area between target onset and 
						the earlier TDT, normalized by that TDT

area_between_tdt_both 	= if both cells select the target, this is the pstch area between the first TDT
						and the second TDT,normalized by that length of time
				
area_after_tdt_xor 		= if only one cell selects the target, this is the pstch area between TDT, 
						and the 90% longest saccade, normalized for that length of time

area_after_tdt_both		= if both cells select the target, this is the pstch area between TDT, 
						and the 90% longest saccade, normalized for that length of time
						
noise_corr				= the noise correlation r value for this pair

pair_nums_all			= the unique pair numbers corresponding to all the pairs in this session
						(NOTE: pair numbers starting with 1 were recorded from F, numbers starting with 
						2 were recorded from Q, and numbers starting with 3 were recorded from S)

pair_nums_no_tdt		= the pair numbers corresponding to pairs in which neither cell selects the target

pair_nums_xor			= the pair numbers corresponding to pairs in which only one cell selects the target
						
pair_nums_both			= the pair numbers corresponding to those pairs in which both cells select the target

fr_before				= the mean firing rate before TDT

fr_after				= the mean firing rate after TDT

peak_fr_no_tdt			= the minimum peak firing rate between cells in a pair where neither selects the target

peak_fr_all				= the minimum peak firing rate between cells in all pairs

rf_title				= the verbose name for this response field (useful for graphs)

hemisphere is coded as follows:
		same hemisphere 		= 1
		different hemisphere	= 0

cell_type is coded as follows:
		vis_vis 			= 1
		vismove_vismove 	= 2
		move_move 			= 3
		vismove_vis 		= 4
		vismove_move		= 5
		vis_move			= 6
		other				= 7