from starsubmit.starsubmit import *
from subprocess import call
from argparse import ArgumentParser


def cli_parser():
    parser = ArgumentParser('cli_parser')
    parser.add_argument('paramater_file', help='File that contains job paramaters')
    parser.add_argument('xml_file', help='Create an xml file with this name')
    parser.add_argument('-s', action='store_true', 
                        help='Submit the xml file after creating')
    return parser

if __name__ == '__main__':
    args = cli_parser().parse_args()

    sub = Submission(args.paramater_file)
    sub.make_xml()
    sub.write_xml(args.xml_file)
    if args.submit:
        call(['star-submit', xml_file])
