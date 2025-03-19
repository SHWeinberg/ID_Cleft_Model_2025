
cluster = parcluster; % get a handle to cluster profile 
cluster.AdditionalProperties.AccountName = 'PAS1622' % set account name 
cluster.AdditionalProperties.WallTime = '00:10:00'; % set wall time to 10 mintues 
cluster.AdditionalProperties.MemUsage =  '6000mb';         
cluster.saveProfile; % locally save the profile



% taskID = getenv('SLURM_ARRAY_TASK_ID')
% feature('numcores')


parpool(cluster,128)

spmd
    fprintf("Worker %d says Hello", spmdIndex);
 end
