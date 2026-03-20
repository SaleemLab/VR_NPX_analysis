function fields_combined = combine_fields_from_multiple_probes(fields1,fields2);
% Functions for combining clusters from multiple probes into one structure
%
all_fields = fieldnames(fields1);


fields_combined = [];

for n = 1:length(all_fields)
    if contains(all_fields{n},'cluster_id')
        for i = 1:length(fields2)
            for ncell = 1:length(fields2(i).cluster_id)
                fields2(i).cluster_id(ncell) = fields2(i).cluster_id(ncell) + 10000;
            end
        end
    end

    for i = 1:length(fields1)
        if contains(all_fields{n},'good_cells')
            fields_combined(i).(all_fields{n}) = [fields1(i).(all_fields{n}) length(fields1(i).cluster_id)+fields2(i).(all_fields{n})];% new index for field2 (i.e. plus field1 cluster numbers)
        elseif  contains(all_fields{n},'shuffled')
            fields_combined(i).(all_fields{n}) = [fields1(i).(all_fields{n}); fields2(i).(all_fields{n})];
        elseif contains(all_fields{n},'x_bin')
            fields_combined(i).(all_fields{n}) = fields1(i).(all_fields{n});
        else
            if size(fields1(i).(all_fields{n}),1)>size(fields1(i).(all_fields{n}),2)
                fields_combined(i).(all_fields{n}) = [fields1(i).(all_fields{n}); fields2(i).(all_fields{n})];
            else
                fields_combined(i).(all_fields{n}) = [fields1(i).(all_fields{n}) fields2(i).(all_fields{n})];
            end
        end
    end
end

end