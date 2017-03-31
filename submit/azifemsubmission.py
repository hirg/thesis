from starsubmit.job import Job
from starsubmit.request import Request
from subprocess import call

do_submit = True

file_list = '$FILELIST'
config_file = 'hbt/femto.config'
macro = 'hbt/AziFemAnalysis.C'
argument_mask = '\\(\\"{}\\",\\"{}\\",{},{},{},{},{},\\"{}\\",\\"{}\\",{}\\)' 
# species = 'UU'
species = 'AuAu'
submission_params_file = 'submit/azifemSubmit.params'
n_events = 99999999

# AuAu - 3341 Files
min_files = [326, 234, 173, 121, 101, 101, 74, 101, 117, 149, 189, 272]
max_files = [335, 239, 176, 140, 120, 108, 99, 108, 124, 152, 197, 279]

# UU - 710 Files
# min_files = [47, 38, 20, 17, 10, 10, 10, 10, 17, 20, 30, 47]
# max_files = [51, 42, 23, 21, 14, 14, 14, 14, 21, 23, 34, 51]

n_zdcBins = 2
n_q2OrMultBins = 2
n_q2MultBins = 5
n_vzBins = 12

for i_zdc in range(0, n_zdcBins):
    for q2OrMult in range(1,n_q2OrMultBins):
        for i_q2MultBins in range(0,n_q2MultBins):
            for i_vz in range(0,n_vzBins):
                if( q2OrMult == 0 ):
                    q2_bin = i_q2MultBins
                    mult_bin = -1
                    bin_label = 'q2_' + str(i_q2MultBins)
                elif( q2OrMult == 1 ):
                    q2_bin = -1
                    mult_bin = i_q2MultBins
                    bin_label = 'mult_' + str(i_q2MultBins)

                job_label = '{}Femto_zdc_{}_{}_Vz_{}'.format(
                                                species, i_zdc, bin_label, i_vz)

                # Set up the command
                # job_label = '%sFemto_zdc_%d_%s_Vz_7' % (species, i_zdc, bin_label)
                out_file = job_label + '_${JOBID}.root'
                log_file = job_label + '_${JOBID}.log'
                xml_file = job_label + '.xml'

                # Hard coded 5cm bin width
                vz_low = -30 + i_vz*5
                vz_high = vz_low + 5

                argument_list = [file_list, out_file, vz_low, vz_high, q2_bin,
                                 mult_bin, i_zdc, config_file, species, n_events]
                argument_string = argument_mask.format(*argument_list)
                command = 'root -b {}{} >& {}'.format(macro, argument_string, log_file)

                # Set up the job/request
                settings = {'minFilesPerProcess':str(min_files[i_vz]), 
                            'maxFilesPerProcess':str(max_files[i_vz])}
                job = Job(submission_params_file)
                job.config['attributes'].update(settings)
                job.add_commands(command)
                request = Request(job)

                # print(request)

                with open(xml_file, 'w') as f:
                    f.write(request.__str__())


                if do_submit:
                    call(['star-submit', xml_file])
