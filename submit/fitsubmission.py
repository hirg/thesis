from starsubmit.starsubmit import *
from subprocess import call


do_submit = True

input_directory = '/star/u/campbell/scratch/200GeV'
config_file = 'fit/fit.config'
fit_macro = 'fit/fitManager.C'
# species = 'UU'
species = 'AuAu'
submission_params_file = 'submit/fitSubmit.params'
n_events = 999999
plus_or_minus = 0

# n_zdcBins = 1
# n_q2OrMultBins = 1
# n_q2MultBins = 1
# n_phiBins = 1
n_zdcBins = 2
n_q2OrMultBins = 2
n_q2MultBins = 5
n_phiBins = 8

fit_options = {'in_file': 'defaultIn.root',
               'out_file': 'defaultOut.root',
               'zdc_bin': 0,
               'plus_minus': 0, # -1 = piMinus, 0 = piCombined, 1 = piPlus
               'q2_or_mult': 0, # 0 = q2, 1 = mult
               'q2Mult_bin': 0,
               'phi_bin': 0,
               'kt_bin': 4, # 4 = Integrated kT
               'fit_range': 0.149,
               'proj_range': 0.03,
               'init_norm': 0.16,
               'init_lambda': 0.45,
               'init_Ro': 27,
               'init_Rs': 22,
               'init_Rl': 32,
               'init_Ros': 0,
               'init_Rol': 0,
               'init_Rsl': 0,
                }

def generate_arg_string(arg_dictionary):
    """ Take a dictionary and returns a a string of arguments """
    ordered_key_list = ['in_file', 'out_file', 'zdc_bin', 'plus_minus', 
        'q2_or_mult', 'q2Mult_bin', 'phi_bin', 'kt_bin', 'fit_range', 
        'proj_range', 'init_norm', 'init_lambda', 'init_Ro', 'init_Rs',
        'init_Rl', 'init_Ros', 'init_Rol', 'init_Rsl']

    arg_string = ''
    for key in ordered_key_list:
        arg_string += str(arg_dictionary[key]) + ','

    return arg_string[:-1] #Ignore the final comma

for i_zdc in range(0, n_zdcBins):
    fit_options['zdc_bin'] = i_zdc

    for q2OrMult in range(0,n_q2OrMultBins):
        fit_options['q2_or_mult'] = q2OrMult

        for i_q2MultBins in range(0,n_q2MultBins):
            fit_options['q2Mult_bin'] = i_q2MultBins
            if( q2OrMult == 0 ):
                q2_bin = i_q2MultBins
                mult_bin = -1
                bin_label = 'q2_' + str(i_q2MultBins)
            elif( q2OrMult == 1 ):
                q2_bin = -1
                mult_bin = i_q2MultBins
                bin_label = 'mult_' + str(i_q2MultBins)

            #Make the xml file
            sub = Submission(submission_params_file)

            # Set up the command
            specifiers = (species, i_zdc, bin_label)
            job_label = '{}Fit_zdc_{}_{}'.format(*specifiers)
            input_label = '{}Femto_zdc_{}_{}'.format(*specifiers)

            fit_options['in_file'] = r'\"{}/{}.root\"'.format(input_directory, input_label)
            fit_options['out_file'] = r'\"{}_{}.root\"'.format(job_label, '${JOBID}')
            log_file = r'{}_{}.log'.format(job_label, '${JOBID}')
            xml_file = job_label + '.xml'



            for i_phiBins in range(0,n_phiBins):
                if i_phiBins:
                    redirect = '>>&'
                else:
                    redirect = '>&'
                fit_options['phi_bin'] = i_phiBins
                args = generate_arg_string(fit_options)
                command = 'root -l -q -b {}\({}\) {} {}'.format(fit_macro, args, redirect, log_file)
                sub.commands.append(command)

            sub.make_xml()
            sub.write_xml(xml_file)

            if do_submit:
                call(['star-submit', xml_file])
