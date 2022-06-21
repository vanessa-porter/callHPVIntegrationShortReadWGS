#!/usr/bin/env python
"""sampletsvtoyaml.py: Converts tab-separated sample information to configuration file for Ploidetect"""

import yaml
import argparse
import os

class inputFileError(ValueError):
    '''raise this when there's a mistake in an input file'''

parser = argparse.ArgumentParser(description = "Converts a tab-separated file of sample information to a samples.yaml file for the Ploidetect pipeline")

parser.add_argument("-t", "--tsv_file", help = "Path to input tsv file", default = "./samples.tsv")
parser.add_argument("-o", "--output_path", help = "Path to config directory. Defaults to ./config/", default = "./config/")
parser.add_argument("-c", "--concat", help = "Should the output be concatenated onto an existing samples.yaml? Accepts True or False")

parser = parser.parse_args()

if not os.path.exists(parser.output_path):
    os.mkdir(parser.output_path)

with open(parser.tsv_file) as f:
    f = [l.strip().split("\t") for l in f.readlines()]

samples = dict(f)
out_dict = {"bams":samples}

if parser.concat == "True":
    out_dict = out_dict["bams"]
    with open(os.path.join(parser.output_path, "samples.yaml"), "r") as f:
        cur_dict = yaml.load(f, Loader = yaml.FullLoader)
        n_dict = {**cur_dict["bams"], **out_dict}
        n_dict = {"bams":n_dict}
    with open(os.path.join(parser.output_path, "samples.yaml"), "w") as f:
        yaml.dump(n_dict, f, default_flow_style=False)
else:
    with open(os.path.join(parser.output_path, "samples.yaml"), "w") as f:
        yaml.dump(out_dict, f, default_flow_style=False)