from starsubmit.job import Job
from starsubmit.request import Request
from subprocess import call
from argparse import ArgumentParser


def cli_parser():
    parser = ArgumentParser('cli_parser')
    parser.add_argument('parameter_file', help='File that contains job parameters')
    parser.add_argument('xml_file', help='Create an xml file with this name')
    parser.add_argument('-s', action='store_true', 
                        help='Submit the xml file after creating')
    return parser

if __name__ == '__main__':
    args = cli_parser().parse_args()

    job = Job(args.parameter_file)
    request = Request(job)

    with open(args.xml_file, 'w') as f:
        f.write(request.__str__())

    if args.s:
        call(['star-submit', args.xml_file])


#     if args.s:
#         call(['star-submit', args.xml_file])
