from starsubmit.starsubmit import *

filelist = '\\"$FILELIST\\"'
n_bins = 5

sub = Submission('200GeVFemto.submit')

# sub.commands = [
#     'echo "Hello World" > hello.log',
#     'ls -lh >> hello.log',
#     'cp hello.log ~/'
#     ]

for i in range(0,n_bins):
    outfile = '"output_${JOBID}.root"'
    command = 'root -l -q -b AziFemAnalysisUU.C\("$FILELIST",' + outfile + ',' + str(i) + '\)'
    sub.commands.append(command)

sub.make_xml()
sub.write_xml()
