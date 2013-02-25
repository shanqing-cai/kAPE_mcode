function raw_fn=getRawFN_(expDir,fn)
[path1,fn1]=fileparts(fn);
[path2,fn2]=fileparts(path1);
[path3,fn3]=fileparts(path2);
[path4,fn4]=fileparts(path3);

raw_fn=fullfile(expDir,fn4,fn3,fn2,fn1);
if ~isequal(raw_fn(end-3:end),'.mat')
    raw_fn=[raw_fn,'.mat'];
end

[ret,hostName]=system('hostname');
hostName=lower(deblank(hostName));

if ~isequal(lower(deblank(hostName)), 'smcg_w510')
    raw_fn = strrep(raw_fn, 'E:', 'D:');
else
    raw_fn = strrep(raw_fn, 'D:', 'E:');
end






return