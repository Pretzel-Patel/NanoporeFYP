'''
pod5_maker.py

This file takes a list of digitised signals and their reads_ids, and stores this information into a POD5 file.

It also shows all the metadata needed to create a POD5 file.
This requires installing the pod5 package, which works on windows.
'''
import uuid
import datetime
import os
import sys
import pod5 as p5
import pytz
import numpy as np
os.chdir(sys.path[0])


# Create 10 random squiggles with random read_ids
signals = []
read_ids = []
num_reads = 10
for i in range(num_reads):
    signals.append(np.random.uniform(200, 800, 4000).astype(np.uint16))
    read_ids.append(f'00000000-0000-0000-0000-{"".join(np.random.randint(0,9,12).astype(str))}')


# Offset for int16 values -  Not important
offset = 10
scale = 0.1755

# Name of new pod5 file
pod5_filename = "new.pod5"
sample_rate = 5000

# Prepare pod5 templates
utc = pytz.timezone('UTC')
pore = p5.Pore(channel=1, well=1, pore_type="not_set")        # I don't think channel and well matter
calibration = p5.Calibration(offset=offset, scale=scale)            # current_pA = scale * (raw_int + offset) = sd_current * current_normalised + mean_current
end_reason = p5.EndReason(reason=p5.EndReasonEnum.UNKNOWN, forced=False)
run_info = p5.RunInfo(acquisition_id='a08e850aaa44c8b56765eee10b386fc3e516a62b', 
                    acquisition_start_time=datetime.datetime(2024, 5, 24, 12, 0, 0,tzinfo=utc), 
                    adc_max=4095, 
                    adc_min=-4096, 
                    context_tags={'basecall_config_filename': 'dna_r9.4.1_450bps_fast.cfg',             # This parameter may need to be changed.
                                    'experiment_duration_set': '180', 
                                    'experiment_type': 'genomic_dna', 
                                    'package': 'bream4', 
                                    'package_version': '4.0.6', 
                                    'sample_frequency': '4000', 
                                    'sequencing_kit': 'sqk-lsk108'}, 
                    experiment_name='', 
                    flow_cell_id='', 
                    flow_cell_product_code='', 
                    protocol_name='c449127e3461a521e0865fe6a88716f6f6b0b30c', 
                    protocol_run_id='df049455-3552-438c-8176-d4a5b1dd9fc5', 
                    protocol_start_time=datetime.datetime(1970, 1, 1, 0, 0,tzinfo=utc), 
                    sample_id='TEST_SAMPLE', 
                    sample_rate=sample_rate, 
                    sequencing_kit='sqk-lsk108', 
                    sequencer_position='MS00000', 
                    sequencer_position_type='minion', 
                    software='python-pod5-converter', 
                    system_name='', 
                    system_type='', 
                    tracking_id={'asic_id': '131070', 
                                'asic_id_eeprom': '0', 
                                'asic_temp': '35.043102', 
                                'asic_version': 'IA02C', 
                                'auto_update': '0', 
                                'auto_update_source': 'https://mirror.oxfordnanoportal.com/software/MinKNOW/', 
                                'bream_is_standard': '0', 
                                'device_id': 'MS00000', 
                                'device_type': 'minion', 
                                'distribution_status': 'modified', 
                                'distribution_version': 'unknown', 
                                'exp_script_name': 'c449127e3461a521e0865fe6a88716f6f6b0b30c', 
                                'exp_script_purpose': 'sequencing_run', 
                                'exp_start_time': '2019-05-13T11:11:43Z', 
                                'flow_cell_id': '', 
                                'guppy_version': '3.0.3+7e7b7d0', 
                                'heatsink_temp': '35.000000', 
                                'hostname': 'happy_fish', 
                                'installation_type': 'prod', 
                                'local_firmware_file': '1', 
                                'operating_system': 'ubuntu 16.04', 
                                'protocol_group_id': 'TEST_EXPERIMENT', 
                                'protocol_run_id': 'df049455-3552-438c-8176-d4a5b1dd9fc5', 
                                'protocols_version': '4.0.6', 
                                'run_id': 'a08e850aaa44c8b56765eee10b386fc3e516a62b', 
                                'sample_id': 'TEST_SAMPLE', 
                                'usb_config': 'MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#ctrl#Auto', 
                                'version': '3.4.0-rc3'})


# POD5 setup - put all reads into a single file
reads = []

for i in range(num_reads):
    # Append read to list of reads
    reads.append(p5.Read(
        read_id=uuid.UUID(read_ids[i]),
        end_reason=end_reason,
        calibration=calibration,
        pore=pore,
        run_info=run_info,
        signal=signals[i],
        read_number=1,
        start_sample=50000000,
        median_before=180.0
    ))

print('Writing reads to pod5 file')
with p5.Writer(pod5_filename) as writer:
    # Write the reads and all of its metadata
    for read in reads:
        writer.add_read(read)

