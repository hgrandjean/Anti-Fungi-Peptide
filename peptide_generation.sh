#!/bin/bash

# Generate a unique run ID
run_id=$(date +"%Y%m%d%H%M%S")

# Generate output file for prints
output_file="$run_id.txt"

# Default verbose mode is false
generate_only=false
selection=true
npep=10
bootstrap=500

# Function to display the help message
show_help() {
    echo "Usage: generate.sh [OPTIONS]"
    echo "Options:"
    echo "  -g             Enable generate-only mode"
    echo "  -s             Enable selection-only mode"
    echo "  --npep x       Set the number of peptides generated (default: 10)"
    echo "  --bootstrap y  Set the number of bootstrap iterations (default: 500)"
    echo "  -h             Show help"
    echo ""
    echo "Example:"
    echo "  generate.sh -v -g --npep 50 --bootstrap 1000"
    echo ""
}

# Parse command-line options
while [[ $# -gt 0 ]]; do
  case $1 in
    -g)
      generate_only=true
      shift
      ;;
    -s)
      selection=false
      shift
      ;;
    --npep)
      npep=$2
      shift 2
      ;;
    --bootstrap)
      bootstrap=$2
      shift 2
      ;;
    -h)
      show_help
      exit 0
      ;;
    *)
      echo "Invalid option: $1" >&2
      show_help
      exit 1
      ;;
  esac
done

# Activate conda environment
# conda activate env

# Shift the positional parameters to skip past the options
shift "$((OPTIND-1))"
if $selection; then
    python generate_peptide.py $npep $bootstrap #> "$output_file" 
  # run Tango on generated file peptide see http://tango.crg.es/about.jsp
  bash results/tango_results/generated_peptides_tango.sh
  mv *_.txt results/tango_results
  grep "AGG" results/tango_results/peptide_agregg.txt | awk '{gsub(/AGG|AMYLO|TURN|HELIX|HELAGG|BETA/,"\t",$0); print;}' > results/tango_results/reformated_aggregated_peptides.tsv
  sed -i '1s/^/\tAGG\tAMYLO\tTURN\tHELIX\tHELAGG\tBETA\n/' results/tango_results/reformated_aggregated_peptides.tsv
  echo ""
  echo "Runned and reformated Tango in vitro prediction"
  echo ""
  if $generate_only ; then 
    exit 1
  fi  
fi


read -e -p "Enter the path to Aggrescan results file: " input_aggrescan_path
# Check if the input file exists
if [ -f "$input_aggrescan_path" ]; then
   python select_active.py  $input_aggrescan_path results/tango_results/reformated_aggregated_peptides.tsv  #> "$output_file"  #check for python file name   
else
  echo "Input file not found: $input_aggrescan_path"
  exit 1
fi
