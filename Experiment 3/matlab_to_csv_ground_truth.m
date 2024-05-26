N = length(verified_reads);

orientations = string(zeros(N,1));
labels = zeros(N,1);
read_ids = string(zeros(N,1));
for i = 1:N
    read = verified_reads(i);
    orientations(i) = 'template';
    if read{1,1}.ground_truth.template == 0
        orientations(i) = 'reverse';
    end
    labels(i) = read{1,1}.ground_truth.tag_idx;
    read_ids(i) = string(read{1,1}.dataset_path(7:42));
end

data = [read_ids, orientations, labels];

writematrix(data, './4_basecalls_csv/adrian_reads.csv')