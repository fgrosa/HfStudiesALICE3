#!/bin/bash

. $ACTS_ROOT/python/setup.sh

export ACTS_SEQUENCER_DISABLE_FPEMON=true

python3 simulate_hf_events.py -c $1