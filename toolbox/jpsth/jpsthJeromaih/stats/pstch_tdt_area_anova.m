function [rf_factor, monkey_factor, data, pair_num, rf_overlap, hemisphere_factor] = pstch_tdt_area_anova(pstch_data, data_type, tdt, timeframe)
	
	if nargin < 4, timeframe = 'all'; end
	
	switch data_type
		case 'actual'
			this_data = pstch_data.actual;
		case 'simulated'
			this_data = pstch_data.sim;
	end

	all = reformat_this_data(this_data.all);
	intersect = reformat_this_data(this_data.int);
	union = reformat_this_data(this_data.union);
	opp = reformat_this_data(this_data.opp);
	xor = reformat_this_data(this_data.xor);
	out = reformat_this_data(this_data.out);
	
	%make all vectors same length (add NaNs where a certain condition doesn't exist) so that we can do a paired wilcoxon test later
    intersect_with_nans = add_missing_vals(all, intersect);
    union_with_nans = add_missing_vals(all, union);
    opp_with_nans = add_missing_vals(all, opp);
    xor_with_nans = add_missing_vals(all, xor);
    out_with_nans = add_missing_vals(all, out);
    
	switch tdt	
		case 'no_tdt'
    		fr_indices = all.peak_fr_no_tdt>=50;
    		all_length = length(all.no_tdt_pstch_area(fr_indices));
    		pair_num = horzcat(all.pair_nums_no_tdt(fr_indices), intersect_with_nans.num_no_tdt(fr_indices), union_with_nans.num_no_tdt(fr_indices), opp_with_nans.num_no_tdt(fr_indices), xor_with_nans.num_no_tdt(fr_indices), out_with_nans.num_no_tdt(fr_indices));
    		data = horzcat(all.no_tdt_pstch_area(fr_indices), intersect_with_nans.no_tdt_area(fr_indices), union_with_nans.no_tdt_area(fr_indices), opp_with_nans.no_tdt_area(fr_indices), xor_with_nans.no_tdt_area(fr_indices), out_with_nans.no_tdt_area(fr_indices));
    		
    	case 'both_tdt'
    		all_length = length(all.area_before_tdt_both);
    		pair_num = horzcat(all.pair_nums_both, intersect_with_nans.num_both,  union_with_nans.num_both, opp_with_nans.num_both, xor_with_nans.num_both, out_with_nans.num_both);
			switch timeframe
				case 'before'
					data = horzcat(all.area_before_tdt_both, intersect_with_nans.before_area_both, union_with_nans.before_area_both, opp_with_nans.before_area_both, xor_with_nans.before_area_both, out_with_nans.before_area_both);
				case 'after'
					data = horzcat(all.area_after_tdt_both, intersect_with_nans.after_area_both, union_with_nans.after_area_both, opp_with_nans.after_area_both, xor_with_nans.after_area_both, out_with_nans.after_area_both);
			end
							
    	case 'xor_tdt'
    		all_length = length(all.area_before_tdt_xor);
    		pair_num = horzcat(all.pair_nums_xor, intersect_with_nans.num_tdt,  union_with_nans.num_tdt, opp_with_nans.num_tdt, xor_with_nans.num_tdt, out_with_nans.num_tdt);
    		switch timeframe
				case 'before'
					data = horzcat(all.area_before_tdt_xor, intersect_with_nans.before_area_xor, union_with_nans.before_area_xor, opp_with_nans.before_area_xor, xor_with_nans.before_area_xor, out_with_nans.before_area_xor);
				case 'after'
					data = horzcat(all.area_after_tdt_xor, intersect_with_nans.after_area_xor, union_with_nans.after_area_xor, opp_with_nans.after_area_xor, xor_with_nans.after_area_xor, out_with_nans.after_area_xor);
            end
	end	
	
	% The monkey factor is coded with Fechner = 1, Quincy = 2, Seymour = 3
	monkey = zeros(1,length(pair_num));
    for i=1:length(pair_num)
       temp = num2str(pair_num(i));
       monkey_factor(i) = str2double(temp(1));
    end
	
	rf_factor = [repmat({'all'},1,all_length), repmat({'intersect'},1,all_length), repmat({'union'},1,all_length), repmat({'opp'},1,all_length), repmat({'xor'},1,all_length), repmat({'out'},1,all_length)];
	rf_factor(isnan(data))=NaN; 
	
	hemisphere_factor = get_hemisphere(pair_num, all.pair_nums_all, all.hemisphere);
	hemisphere_factor(isnan(data))=NaN;
	
	rf_overlap = overlap_nonoverlap_factor(pair_num, rf_factor);
	rf_overlap(isnan(data))=NaN;
end

function [with_nans] = add_missing_vals(all,rf)
	for i=1:length(all.pair_nums_xor)
		if all.pair_nums_xor(i) ~= rf.pair_nums_xor
			this_num_tdt(i) = NaN;
			this_before_area(i) = NaN;
			this_after_area(i) = NaN;
		else
			this_num_tdt(i) = all.pair_nums_xor(i);
			this_before_area(i) = rf.area_before_tdt_xor(rf.pair_nums_xor==all.pair_nums_xor(i));
			this_after_area(i) = rf.area_after_tdt_xor(rf.pair_nums_xor==all.pair_nums_xor(i));
		end
	end
	
	for i=1:length(all.pair_nums_both)
		if all.pair_nums_both(i) ~= rf.pair_nums_both
			this_num_both(i) = NaN;
			this_before_area_both(i) = NaN;
			this_after_area_both(i) = NaN;
		else
			this_num_both(i) = all.pair_nums_both(i);
			this_before_area_both(i) = rf.area_before_tdt_both(rf.pair_nums_both==all.pair_nums_both(i));
			this_after_area_both(i) = rf.area_after_tdt_both(rf.pair_nums_both==all.pair_nums_both(i));
		end
	end
	
	for i=1:length(all.pair_nums_no_tdt)
		if all.pair_nums_no_tdt(i) ~= rf.pair_nums_no_tdt
			this_num_no_tdt(i) = NaN;
			this_area_no_tdt(i) = NaN;
		else
			this_num_no_tdt(i) = all.pair_nums_no_tdt(i);
			this_area_no_tdt(i) = rf.no_tdt_pstch_area(rf.pair_nums_no_tdt==all.pair_nums_no_tdt(i));
		end
	end
	
	with_nans = struct('num_tdt',this_num_tdt,'num_both',this_num_both,'before_area_xor',this_before_area,...
						'after_area_xor',this_after_area,'before_area_both',this_before_area_both,'after_area_both',this_after_area_both,...
						'no_tdt_area',this_area_no_tdt,'num_no_tdt', this_num_no_tdt);
	
end

function [rf_overlap] = overlap_nonoverlap_factor(pair_num, rf_factor)

	for i=1:length(pair_num)
		rf = rf_factor(i);
		if strcmp(rf,'intersect')
			pair_num_intersect(i) = pair_num(i);
		end 
    end
	
    pair_num_intersect = pair_num_intersect(~isnan(pair_num_intersect)& pair_num_intersect>=1);
    
    for i=1:length(pair_num)
        if intersect(pair_num(i),pair_num_intersect)
            rf_overlap(i) = {'overlap'};
        elseif isnan(pair_num(i))
            rf_overlap(i) = {'NaN'};
        else                
            rf_overlap(i) = {'nonoverlap'};
        end
    end
end

function [hemisphere_factor] = get_hemisphere(pair_num, pair_nums_all, hemisphere)
	hemisphere_factor = zeros(1,length(pair_num));
	for i=1:length(pair_num)
		pair_index = find(pair_nums_all==pair_num(i));
		if ~isempty(pair_index)
			hemisphere_factor(i) = hemisphere(pair_nums_all==pair_num(i));
		else
			hemisphere_factor(i) = NaN;
		end
    end
end

function [formatted_data] = reformat_this_data(this_data)

	formatted_data.no_tdt_pstch_area = []; formatted_data.area_before_tdt_xor = []; formatted_data.area_before_tdt_both = [];
	formatted_data.area_after_tdt_xor = []; formatted_data.area_after_tdt_both = []; formatted_data.pair_nums_no_tdt = []; 
	formatted_data.pair_nums_xor = []; formatted_data.pair_nums_both = []; formatted_data.peak_fr_no_tdt = []; 
	formatted_data.hemisphere = []; formatted_data.pair_nums_all = [];
    
	for i=1:length(this_data)
		formatted_data.no_tdt_pstch_area = [formatted_data.no_tdt_pstch_area, this_data(i).no_tdt_pstch_area(~isnan(this_data(i).no_tdt_pstch_area))];
		formatted_data.area_before_tdt_xor = [formatted_data.area_before_tdt_xor, this_data(i).area_before_tdt_xor(~isnan(this_data(i).area_before_tdt_xor))];
		formatted_data.area_before_tdt_both = [formatted_data.area_before_tdt_both, this_data(i).area_before_tdt_both(~isnan(this_data(i).area_before_tdt_both))];
		formatted_data.area_after_tdt_xor = [formatted_data.area_after_tdt_xor, this_data(i).area_after_tdt_xor(~isnan(this_data(i).area_after_tdt_xor))];
		formatted_data.area_after_tdt_both = [formatted_data.area_after_tdt_both, this_data(i).area_after_tdt_both(~isnan(this_data(i).area_after_tdt_both))];
		formatted_data.pair_nums_no_tdt = [formatted_data.pair_nums_no_tdt, this_data(i).pair_nums_no_tdt(~isnan(this_data(i).pair_nums_no_tdt))];
		formatted_data.pair_nums_xor = [formatted_data.pair_nums_xor, this_data(i).pair_nums_xor(~isnan(this_data(i).pair_nums_xor))];
		formatted_data.pair_nums_both = [formatted_data.pair_nums_both, this_data(i).pair_nums_both(~isnan(this_data(i).pair_nums_both))];
		formatted_data.peak_fr_no_tdt = [formatted_data.peak_fr_no_tdt, this_data(i).peak_fr_no_tdt(~isnan(this_data(i).peak_fr_no_tdt))];
		formatted_data.hemisphere = [formatted_data.hemisphere, this_data(i).hemisphere(~isnan(this_data(i).hemisphere))];
		formatted_data.pair_nums_all = [formatted_data.pair_nums_all, this_data(i).pair_nums_all(~isnan(this_data(i).pair_nums_all))];
	end
	
end