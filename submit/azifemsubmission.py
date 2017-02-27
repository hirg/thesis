from starsubmit.starsubmit import *
from subprocess import call


do_submit = True

file_list = '$FILELIST'
config_file = 'hbt/femto.config'
macro = 'hbt/AziFemAnalysis.C'
# species = 'UU'
species = 'AuAu'
submission_params_file = 'submit/azifemSubmit.params'
n_events = 99999999

n_zdcBins = 1
n_q2OrMultBins = 2
n_q2MultBins = 5

for i_zdc in range(0, n_zdcBins):
    for q2OrMult in range(1,n_q2OrMultBins):
        for i_q2MultBins in range(0,n_q2MultBins):
            if( q2OrMult == 0 ):
                q2_bin = i_q2MultBins
                mult_bin = -1
                bin_label = 'q2_' + str(i_q2MultBins)
            elif( q2OrMult == 1 ):
                q2_bin = -1
                mult_bin = i_q2MultBins
                bin_label = 'mult_' + str(i_q2MultBins)

            # Set up the command
            job_label = '%sFemto_zdc_%d_%s_Vz_7' % (species, i_zdc, bin_label)
            out_file = job_label + '_${JOBID}.root'
            log_file = job_label + '_${JOBID}.log'
            xml_file = job_label + '.xml'

            argsTuple = (file_list, out_file, q2_bin, mult_bin, i_zdc,
                            config_file, species, n_events)
            args = ( '\\"%s\\",\\"%s\\",%d,%d,%d,\\"%s\\",\\"%s\\",%d'
                        % argsTuple )
            command = 'root -l -q -b %s\(%s\) >& %s' % (macro, args, log_file)

            #Make the xml file
            sub = Submission(submission_params_file)
            sub.commands.append(command)

            sub.make_xml()
            sub.write_xml(xml_file)

            if do_submit:
                call(['star-submit', xml_file])
