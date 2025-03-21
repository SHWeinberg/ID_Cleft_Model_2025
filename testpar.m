
walltime = getenv('EXP_WALLTIME');

cluster = parcluster; % get a handle to cluster profile 
cluster.AdditionalProperties.AccountName = 'PAS1622' % set account name 
cluster.AdditionalProperties.WallTime = walltime; % set wall time to 10 mintues 
cluster.AdditionalProperties.MemPerCPU = '10gb'; 
cluster.saveProfile; % locally save the profile



% taskID = getenv('SLURM_ARRAY_TASK_ID')
% feature('numcores')


parpool(cluster,2)

spmd

    fprintf("Worker %d says Hello", spmdIndex);

    localDir = getenv('TMPDIR') + "/";
    disp(localDir)
    scratchDir = '/fs/scratch/PAS1622/nickmoise/ID_2025/';

    save_name_data = localDir + string(spmdIndex) + ".mat";
    scratch_save_name = scratchDir + string(spmdIndex)+ ".mat";

    save_data_test(save_name_data, spmdIndex)
    copyfile(save_name_data, scratch_save_name);

end



function save_data_test(save_name_data,spmdIndex)
    save(save_name_data,'spmdIndex');
end
