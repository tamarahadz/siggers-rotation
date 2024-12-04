#!/bin/bash

# Some common functionality to avoid repeating it in each qsub

function prerun() {
  echo "=========================================================="
  echo "Starting on       : $(date)"
  echo "Run by            : $(whoami)"
  echo "Running on node   : $(hostname)"
  echo "Current directory : $(pwd)"
  echo "Current job ID    : $JOB_ID"
  echo "Current job name  : $JOB_NAME"
  echo "Task index number : $SGE_TASK_ID"
  echo "=========================================================="
}

function run() {
  echo "Running command: '$1'"
  eval $1
}

function postrun() {
  echo "=========================================================="
  echo "Finished on       : $(date)"
  echo "=========================================================="
}