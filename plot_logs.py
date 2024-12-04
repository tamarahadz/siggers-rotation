#!/usr/bin/env python3

from argparse import ArgumentParser
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import sys

def plot_logs(infile: Path, outfile: Path):

    logfile = pd.read_csv(infile, sep='\t')
    pearsons = logfile[["Iteration", "Validation Count Pearson","Validation Profile Pearson"]]

    plt.plot(pearsons["Iteration"], pearsons["Validation Count Pearson"], label='Count')
    plt.plot(pearsons["Iteration"], pearsons["Validation Profile Pearson"], label='Profile')
    plt.title("Pearson Correlation During Validation")
    plt.legend()

    plt.savefig(outfile)

if __name__ == '__main__':
  argparser = ArgumentParser()
  argparser.add_argument("-i", dest="infile", type=Path)
  argparser.add_argument("-o", dest="outfile", type=Path)
  args = argparser.parse_args()

  plot_logs(args.infile, args.outfile)
