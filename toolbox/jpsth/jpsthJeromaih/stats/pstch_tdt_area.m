function [pstch_data] = pstch_tdt_area(directory_path)
	
	[files] = dir(fullfile(directory_path,'jpsth_*.mat'));

  	for i=1:length(files)
		load(fullfile(directory_path,files(i).name))
		if strcmp(files(i).name(7),'f')
			pair_num_prefix = '1';
		elseif strcmp(files(i).name(7),'q')
			pair_num_prefix = '2';
		elseif strcmp(files(i).name(7),'s')
			pair_num_prefix = '3';
		end
	
        pstch_data.actual.all(i) = extract_data(data,'Target in All Loca','actual', pair_num_prefix, i);
		pstch_data.actual.int(i) = extract_data(data,'Target In Intersec','actual', pair_num_prefix, i);
		pstch_data.actual.union(i) = extract_data(data, 'Target In United R','actual', pair_num_prefix, i);
		pstch_data.actual.opp(i) = extract_data(data, 'Target In Opposite','actual', pair_num_prefix, i);
		pstch_data.actual.xor(i) = extract_data(data, 'Target In One Neur','actual', pair_num_prefix, i);
		pstch_data.actual.out(i) = extract_data(data, 'Target Outside Res','actual', pair_num_prefix, i);
		
		pstch_data.sim.all(i) = extract_data(data,'Target in All Loca','simulated', pair_num_prefix, i);
		pstch_data.sim.int(i) = extract_data(data,'Target In Intersec','simulated', pair_num_prefix, i);
		pstch_data.sim.union(i) = extract_data(data, 'Target In United R','simulated', pair_num_prefix, i);
		pstch_data.sim.opp(i) = extract_data(data, 'Target In Opposite','simulated', pair_num_prefix, i);
		pstch_data.sim.xor(i) = extract_data(data, 'Target In One Neur','simulated', pair_num_prefix, i);
		pstch_data.sim.out(i) = extract_data(data, 'Target Outside Res','simulated', pair_num_prefix, i);
		
	end
  
end

function [rf_struct] = extract_data(data, rf_name, data_type, pair_num_prefix, i)
	
    pair = 0;
	tdt_pair = 0;
	both_tdt_pair = 0;
	no_tdt_pair = 0;
	noise_corr = [];
	pair_nums_all = [];
    pair_nums_no_tdt = [];
    pair_nums_xor = [];
    pair_nums_both = [];
    signal_list = [];
    pstch_area_all = NaN;
    pstch_vector = [];
    no_tdt_pstch_area = NaN;
    area_before_tdt = NaN;
    area_after_tdt = NaN;
    area_before_tdt_both = NaN;
    area_between_tdt_both = NaN;
	area_after_tdt_both = NaN;
    fr_before = NaN;
    fr_after = NaN;
    pair_min_peak_fr_no_tdt = NaN;
    pair_min_peak_fr_all = NaN;
    hemisphere = [];
    cell_type = [];
    
    [this_pair_array] = find_pairs(data);

	for j=1:length(data.pair)
			rf = data.pair(j).field.name;
			time_window = data.alignments(data.pair(j).align).time_window;
			binwidth = data.binwidth;
			time_vector = time_window(1)+(binwidth/2):binwidth:time_window(2)-(binwidth/2);
			
		%select only comparisons using the specified response field and
		%select only Target aligned pairs (Saccade aligned does not make sense in comparison involving TDT)
		if strcmp(rf(1:18),rf_name) && data.pair(j).align == 1 
			pair = pair+1;
			
			%same hemisphere coded as 1, different hemispheres coded as 0
			if strcmp(data.signals(data.pair(j).s1).Hemisphere(1),data.signals(data.pair(j).s2).Hemisphere(1))
				hemisphere(pair) = 1;
			else
				hemisphere(pair) = 0;
			end
			
			% cell type judgements based on numerical ratings (Visual and Move from spreadsheet)
			% vis_vis 			= 1
			% vismove_vismove 	= 2
			% move_move 		= 3
			% vismove_vis 		= 4
			% vismove_move		= 5
			% vis_move			= 6
			% other				= 7
			if str2double(char(data.signals(data.pair(j).s1).Visual)) >= 0.5 && str2double(char(data.signals(data.pair(j).s1).Move)) < 0.5 ...
				&& str2double(char(data.signals(data.pair(j).s2).Visual)) >= 0.5 && str2double(char(data.signals(data.pair(j).s2).Move)) < 0.5
					cell_type(pair) = 1;
			elseif str2double(char(data.signals(data.pair(j).s1).Visual)) >= 0.5 && str2double(char(data.signals(data.pair(j).s1).Move)) >= 0.5 ...
				&& str2double(char(data.signals(data.pair(j).s2).Visual)) >= 0.5 && str2double(char(data.signals(data.pair(j).s2).Move)) >= 0.5
					cell_type(pair) = 2;
			elseif str2double(char(data.signals(data.pair(j).s1).Visual)) < 0.5 && str2double(char(data.signals(data.pair(j).s1).Move)) >= 0.5 ...
				&& str2double(char(data.signals(data.pair(j).s2).Visual)) < 0.5 && str2double(char(data.signals(data.pair(j).s2).Move)) >= 0.5
					cell_type(pair) = 3;
			elseif (str2double(char(data.signals(data.pair(j).s1).Visual)) >= 0.5 && str2double(char(data.signals(data.pair(j).s1).Move)) >= 0.5 ...
				&& str2double(char(data.signals(data.pair(j).s2).Visual)) >= 0.5 && str2double(char(data.signals(data.pair(j).s2).Move)) < 0.5) ...
				|| (str2double(char(data.signals(data.pair(j).s1).Visual)) >= 0.5 && str2double(char(data.signals(data.pair(j).s1).Move)) < 0.5 ...
				&& str2double(char(data.signals(data.pair(j).s2).Visual)) >= 0.5 && str2double(char(data.signals(data.pair(j).s2).Move)) >= 0.5)
					cell_type(pair) = 4;
			elseif (str2double(char(data.signals(data.pair(j).s1).Visual)) >= 0.5 && str2double(char(data.signals(data.pair(j).s1).Move)) >= 0.5 ...
				&& str2double(char(data.signals(data.pair(j).s2).Visual)) < 0.5 && str2double(char(data.signals(data.pair(j).s2).Move)) >= 0.5) ...
				|| (str2double(char(data.signals(data.pair(j).s1).Visual)) < 0.5 && str2double(char(data.signals(data.pair(j).s1).Move)) >= 0.5 ...
				&& str2double(char(data.signals(data.pair(j).s2).Visual)) >= 0.5 && str2double(char(data.signals(data.pair(j).s2).Move)) >= 0.5)
					cell_type(pair) = 5;
			elseif (str2double(char(data.signals(data.pair(j).s1).Visual)) >= 0.5 && str2double(char(data.signals(data.pair(j).s1).Move)) < 0.5 ...
				&& str2double(char(data.signals(data.pair(j).s2).Visual)) < 0.5 && str2double(char(data.signals(data.pair(j).s2).Move)) >= 0.5) ...
				|| (str2double(char(data.signals(data.pair(j).s1).Visual)) < 0.5 && str2double(char(data.signals(data.pair(j).s1).Move)) >= 0.5 ...
				&& str2double(char(data.signals(data.pair(j).s2).Visual)) >= 0.5 && str2double(char(data.signals(data.pair(j).s2).Move)) < 0.5)
					cell_type(pair) = 6;
			else
					cell_type(pair) = 7;
			end
			
			noise_corr(pair) = data.result(j).actual.noise_corr.noise_corr_r;
			pair_min_peak_fr_all(pair) = min(max(data.result(j).actual.psth_1/(binwidth/1000)),max(data.result(j).actual.psth_2/(binwidth/1000)));
			
			switch data_type
				case 'actual'
					pstch_vector(pair,:) = [data.result(j).actual.pstch, NaN*ones(1,1200-length(data.result(j).actual.pstch))];
					pstch_area_all(pair) = binwidth*sum(data.result(j).actual.pstch(find(time_vector>=0,1):end))/length(data.result(j).actual.pstch(find(time_vector>=0,1):end));
				case 'simulated'
					pstch_vector(pair,:) = [data.result(j).simulated.pstch, NaN*ones(1,1200-length(data.result(j).simulated.pstch))];
					pstch_area_all(pair) = binwidth*sum(data.result(j).simulated.pstch(find(time_vector>=0,1):end))/length(data.result(j).simulated.pstch(find(time_vector>=0,1):end));
			end
			
			%finding a unique number to apply to each pair
			for m=1:size(this_pair_array,1)
				if strcmp(this_pair_array(m,1), data.signals(data.pair(j).s1).name) && strcmp(this_pair_array(m,2), data.signals(data.pair(j).s2).name)
					pair_nums_all(pair) = str2double(strcat(pair_num_prefix,num2str(i),num2str(m)));
				end
			end
			
			tdt1 = str2double(cell2mat(data.signals(data.pair(j).s1).TDT));
			tdt2 = str2double(cell2mat(data.signals(data.pair(j).s2).TDT));
						
			if isempty(tdt1), tdt1 = NaN; end
			if isempty(tdt2), tdt2 = NaN; end
			
			if (~isnan(tdt1) && isnan(tdt2)) || (~isnan(tdt2) && isnan(tdt1))
				tdt_pair = tdt_pair + 1;
				tdt = min(tdt1, tdt2);
				tdt_index = find(time_vector>=floor(tdt),1);
				
				[signal_list, firing_rate1_before_tdt(tdt_pair), firing_rate1_after_tdt(tdt_pair)] = extract_spike_rate(data, j, signal_list, tdt1, '1', binwidth, time_vector, tdt_index);
				[signal_list, firing_rate2_before_tdt(tdt_pair), firing_rate2_after_tdt(tdt_pair)] = extract_spike_rate(data, j, signal_list, tdt2, '2', binwidth, time_vector, tdt_index);
				
				fr_before = [firing_rate1_before_tdt(~isnan(firing_rate1_before_tdt)), firing_rate2_before_tdt(~isnan(firing_rate2_before_tdt))];
				fr_after = [firing_rate1_after_tdt(~isnan(firing_rate1_after_tdt)), firing_rate2_after_tdt(~isnan(firing_rate2_after_tdt))];

				pair_nums_xor(tdt_pair) = pair_nums_all(pair);
				
				switch data_type
					case 'actual'
						area_before_tdt(tdt_pair) = binwidth*sum(data.result(j).actual.pstch(find(time_vector>=0,1):tdt_index-1))/length(data.result(j).actual.pstch(find(time_vector>=0,1):tdt_index-1));
						area_after_tdt(tdt_pair) = binwidth*sum(data.result(j).actual.pstch(tdt_index:end))/length(data.result(j).actual.pstch(tdt_index:end));
					case 'simulated'
						area_before_tdt(tdt_pair) = binwidth*sum(data.result(j).simulated.pstch(find(time_vector>=0,1):tdt_index-1))/length(data.result(j).simulated.pstch(find(time_vector>=0,1):tdt_index-1));
						area_after_tdt(tdt_pair) = binwidth*sum(data.result(j).simulated.pstch(tdt_index:end))/length(data.result(j).simulated.pstch(tdt_index:end));
				end
				
			elseif ~isnan(tdt1) && ~isnan(tdt2)
				both_tdt_pair = both_tdt_pair +1;
				pair_nums_both(both_tdt_pair) = pair_nums_all(pair);
				tdt_index1 = find(time_vector>=floor(tdt1),1);
				tdt_index2 = find(time_vector>=floor(tdt2),1);
				
				switch data_type
					case 'actual'
						area_before_tdt_both(both_tdt_pair) = binwidth*sum(data.result(j).actual.pstch(find(time_vector>=0,1):tdt_index1-1))/length(data.result(j).actual.pstch(find(time_vector>=0,1):tdt_index1-1));
						area_between_tdt_both(both_tdt_pair) = binwidth*sum(data.result(j).actual.pstch(tdt_index1:tdt_index2-1))/length(data.result(j).actual.pstch(tdt_index1:tdt_index2-1));
						area_after_tdt_both(both_tdt_pair) = binwidth*sum(data.result(j).actual.pstch(tdt_index2:end))/length(data.result(j).actual.pstch(tdt_index2:end));
					case 'simulated'
						area_before_tdt_both(both_tdt_pair) = binwidth*sum(data.result(j).simulated.pstch(find(time_vector>=0,1):tdt_index1-1))/length(data.result(j).simulated.pstch(find(time_vector>=0,1):tdt_index1-1));
						area_between_tdt_both(both_tdt_pair) = binwidth*sum(data.result(j).simulated.pstch(tdt_index1:tdt_index2-1))/length(data.result(j).simulated.pstch(tdt_index1:tdt_index2-1));
						area_after_tdt_both(both_tdt_pair) = binwidth*sum(data.result(j).simulated.pstch(tdt_index2:end))/length(data.result(j).simulated.pstch(tdt_index2:end));
				end
			
			
			else
				no_tdt_pair = no_tdt_pair + 1;
				pair_nums_no_tdt(no_tdt_pair) = pair_nums_all(pair);
				pair_min_peak_fr_no_tdt(no_tdt_pair) = min(max(data.result(j).actual.psth_1/(binwidth/1000)),max(data.result(j).actual.psth_2/(binwidth/1000)));
					
				switch data_type
					case 'actual'
						no_tdt_pstch_area(no_tdt_pair) = binwidth*sum(data.result(j).actual.pstch(find(time_vector>=0,1):end))/length(data.result(j).actual.pstch(1:end));
					case 'simulated'
						no_tdt_pstch_area(no_tdt_pair) = binwidth*sum(data.result(j).simulated.pstch(find(time_vector>=0,1):end))/length(data.result(j).simulated.pstch(1:end));
				end
				
			end	
		end
	end
	
	switch rf_name
		case 'Target in All Loca'
			rf_title = 'Target in All RF';
		case 'Target In Intersec'
			rf_title = 'Target in Pair RF (Intersection)';
		case 'Target In United R'
			rf_title = 'Target in Pair RF (Union)';
		case 'Target In Opposite'
			rf_title = 'Distractor in Pair RF';
		case 'Target In One Neur'
			rf_title = 'Target In One Neuron''s Response Field of Only';
		case 'Target Outside Res'
			rf_title = 'Target Outside Response Field of Both Neurons';
	end	

	%put the data into a struct so it can be accessed later by other functions
	rf_struct = struct('pstch_vector',pstch_vector,'pstch_area_all',pstch_area_all,'no_tdt_pstch_area',no_tdt_pstch_area, 'area_before_tdt_xor', area_before_tdt, ...
                        'area_before_tdt_both', area_before_tdt_both, 'area_between_tdt_both', area_between_tdt_both, 'area_after_tdt_xor', area_after_tdt, 'area_after_tdt_both',area_after_tdt_both,...
                        'noise_corr', noise_corr, 'pair_nums_all', pair_nums_all,'pair_nums_no_tdt', pair_nums_no_tdt, 'pair_nums_xor', pair_nums_xor, ...
                        'pair_nums_both',pair_nums_both,'fr_before', fr_before,'fr_after', fr_after, 'peak_fr_no_tdt', pair_min_peak_fr_no_tdt, 'peak_fr_all', ...
                        pair_min_peak_fr_all, 'rf_title', rf_title, 'hemisphere',hemisphere, 'cell_type', cell_type);
end

function [signal_list, firing_rate_before_tdt, firing_rate_after_tdt] = extract_spike_rate(data, j, signal_list, tdt, signal_num, binwidth, time_vector, tdt_index)
    firing_rate_before_tdt = NaN;
    firing_rate_after_tdt = NaN;  
    
	switch signal_num
		case '1'
			name = data.signals(data.pair(j).s1).name;
			psth = data.result(j).actual.psth_1;
		case '2'
			name = data.signals(data.pair(j).s2).name;
			psth = data.result(j).actual.psth_2;
			
	end
	
	[unique_signal_num, signal_list] = unique_signal_check(signal_list, name);
	
	if ~isempty(tdt) && ~isnan(tdt)
		if unique_signal_num	
			firing_rate_before_tdt = mean(psth(find(time_vector>=0,1):tdt_index-1))/(binwidth/1000);
			firing_rate_after_tdt = mean(psth(tdt_index:end))/(binwidth/1000);
		end
	end
end