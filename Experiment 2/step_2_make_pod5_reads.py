'''
Experiment 1 Step 2
step_2_make_pod5_reads.py
'''
import uuid
import csv
import datetime
from itertools import product
import os
import sys
os.chdir(sys.path[0])

from step_0_configure import num_trials, sd_min, sd_max, dwell_min, dwell_max, choice_list, offset, scale, model, dorado_call_auto, sample_rate, squiggle_generation

# Get saved data
def load_list_of_lists_from_csv(filename):
    with open(filename, 'r') as csvfile:
        csv_reader = csv.reader(csvfile)
        data = [list(map(int, row)) for row in csv_reader]
    return data


if squiggle_generation == 'scrappie':
    print('No pod5 file needs to be generated since scrappie python package can be used.\n')
    if dorado_call_auto:
        import step_3_summarise
    sys.exit()
else:
    #POD5 maker
    import pod5 as p5
    import pytz

# Check if file exists
pod5_filename = "4_reads/all_reads.pod5"
command = f'dorado basecaller -v --emit-fastq --batchsize 64 ../models_dorado/{model} 4_reads/all_reads.pod5 > 5_rx_msg/all_calls.fa'
if os.path.exists(pod5_filename):
    delete = int(input('all_reads.pod5 exists. Would you like to delete? Enter (0) or (1): '))
    if delete:
        # Delete pod5 file and continue pod5 creation process
        os.remove(pod5_filename)
        print('all_reads.pod5 has been deleted.')
    else:
        # Keep pod5 file and cancel pod5 creation process, instead providing the required command.
        print(f'Run this in windows shell\n{command}')
        sys.exit()

# Prepare pod5 templates
utc = pytz.timezone('UTC')
pore = p5.Pore(channel=1, well=1, pore_type="not_set")        # I don't think channel and well matter
calibration = p5.Calibration(offset=offset, scale=scale)            # current_pA = scale * (raw_int + offset) = sd_current * current_normalised + mean_current
end_reason = p5.EndReason(reason=p5.EndReasonEnum.UNKNOWN, forced=False)
run_info = p5.RunInfo(acquisition_id='a08e850aaa44c8b56765eee10b386fc3e516a62b', 
                    acquisition_start_time=datetime.datetime(2019, 5, 13, 11, 11, 43,tzinfo=utc), 
                    adc_max=4095, 
                    adc_min=-4096, 
                    context_tags={'basecall_config_filename': 'dna_r9.4.1_450bps_fast.cfg', 
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
    # unsure about adc_max and adc_min


# POD5 setup - put all reads into a single file
reads = []

print('cb | sd | dwell')
for choice, sd_mean_scaled, dwell_mean in product(choice_list, range(sd_min,sd_max+1), range(dwell_min,dwell_max+1)):
    
    # File names
    output_file = f'2_tx_msg/tx_message_sd{sd_mean_scaled}_K{dwell_mean}_cb{choice}.fa'
    csv_filename = f'3_squiggles/squiggles_sd{sd_mean_scaled}_K{dwell_mean}_cb{choice}.csv'
    
    print(f'{choice:2.0f}   {sd_mean_scaled:2.0f}   {dwell_mean:2.0f}')

    sd_mean = sd_mean_scaled / 100

    all_signals = load_list_of_lists_from_csv(csv_filename)
    assert num_trials == len(all_signals)

    for i in range(num_trials):
        signal_int16 = all_signals[i]
        read_id = uuid.UUID(f"{str(sd_mean_scaled).zfill(8)}-{str(dwell_mean).zfill(4)}-{str(choice).zfill(4)}-0000-0000{str(i).zfill(8)}")
        # Process 5a. POD5 maker
        reads.append(p5.Read(
            read_id=read_id,
            end_reason=end_reason,
            calibration=calibration,
            pore=pore,
            run_info=run_info,
            signal=signal_int16,
            read_number=1,
            start_sample=50000000,
            median_before=180.0
        ))

print('Writing reads to pod5 file')
with p5.Writer(pod5_filename) as writer:
    # Write the reads and all of its metadata
    for read in reads:
        writer.add_read(read)

# Run dorado in command line using python
if dorado_call_auto:
    os.system(command)
    import step_3_summarise
else:
    print(f'Run this in windows shell\n{command}')
