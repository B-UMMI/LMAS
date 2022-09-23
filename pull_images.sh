#!/bin/bash

if [ -x "$(command -v shifter)" ]; then
  echo "Pulling images with shifter..." >&2
  pull_command="shifterimg pull "

elif [ -x "$(command -v singularity)" ]; then
  echo "Pulling images with singularity..." >&2
  pull_command="singularity pull docker://"

elif [ -x "$(command -v docker)" ]; then
  echo "Pulling images with docker..." >&2
  pull_command="docker pull "

else
  echo "No container software found. Exiting..." >&2
fi

eval $pull_command"pcerqueira/bbtools:38.44"
eval $pull_command"cimendes/abyss:2.3.1-1"
eval $pull_command"cimendes/gatb-minia-pipeline:31.07.2020-1"
eval $pull_command"cimendes/idba:1.1.3-1"
eval $pull_command"cimendes/megahit-assembler:1.2.9-1"
eval $pull_command"cimendes/mhm2:v2.0.0-65-gaad446d-generic"
eval $pull_command"cimendes/spades:3.15.3-1"
eval $pull_command"cimendes/minia:3.2.6-1"
eval $pull_command"cimendes/skesa:2.5.0-1"
eval $pull_command"cimendes/unicycler:0.4.9-1"
eval $pull_command"cimendes/velvetoptimiser:2.2.6-1"
eval $pull_command"cimendes/minimap2:2.22-1"
