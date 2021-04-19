%%
ppn = 32; 
totalProcs = 31;
ClusterInfo.setQueueName('matlab');
ClusterInfo.setProcsPerNode(ppn);
c = parcluster;
%%

% when you want to run a job using 2K digraphs and 30K parameters
[noneighbor1] = c.batch(@modelAll_v3_maptorange,1,{randDigraphs2000, newParameters2},'Pool',totalProcs,'AttachedFiles',{'/Users/aylab/Documents/MATLAB/SegmentationClockCode_Summer2020/SegClockNoNeighbor'});

