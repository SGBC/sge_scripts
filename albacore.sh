#!/usr/bin/env bash

# basecalling for nanopore reads (with albacore)

usage="$(basename "$0") [-h] [-e email] [-c cpus] [-m ram] [-t time]
    [-i fast5 dir] [-o output] [-f flowcell] [-k kit] --

        basecalling for nanopore reads (with albacore)

        qsub options:
    		-h show this useful help
            -e <string> your email address
            -c <int> number of cpus to use
            -m <int> number of ram per cpu to use
            -t <time> upper time limit for the job, in hours

        albacore options:
    		-i <dir> directory containing read fast5 files
            -o <path> path to save output
            -f <string> flowcell used during the sequencing run
            -k <string> kit used during the sequencing run"

while getopts "he:c:m:t:i:s:f:k:o:" option
    do
        case "$option" in
            h) echo "$usage"
            exit 0
                ;;
            e) email="$OPTARG"
                ;;
            c) cpus="$OPTARG"
                ;;
            m) mem="$OPTARG"
                ;;
            t) run_time="$OPTARG"
                ;;
            i) directory=$(readlink -f "$OPTARG")
                ;;
            o) output=$(readlink -f "$OPTARG")
                ;;
            f) flowcell="$OPTARG"
                ;;
            k) kit="$OPTARG"
                ;;
            :) printf "missing argument for -%s\n" "$OPTARG" >&2
            echo "$usage" >&2
            exit 1
                ;;
            \?) printf "illegal option -%s\n" "$OPTARG" >&2
            echo "$usage">&2
            exit 1
                ;;
        esac
    done
shift $((OPTIND-1))

# Checking if the user is executing the script on the server
if ! [[ $(uname -n) =~ "planetsmasher" ]]
    then
        printf "WARNING: You are not connected to planetsmasher!\n \
        Connect to the server and try again.\n" >&2
        exit 1
fi

echo "Checking parameters"
if [ -z "${email+x}" ]
    then printf "%s\n\nEmply value for -e\n" "$usage" >&2; exit 1
fi
if [ -z "${cpus+x}" ]
    then printf "%s\n\nEmply value for -c\n" "$usage" >&2; exit 1
fi
if [ -z "${mem+x}" ]
    then printf "%s\n\nEmply value for -m\n" "$usage" >&2; exit 1
fi
if [ -z "${run_time+x}" ]
    then printf "%s\n\nEmply value for -t\n" "$usage" >&2; exit 1
fi
if [ -z "${directory+x}" ]
    then printf "%s\n\nEmply value for -i\n" "$usage" >&2; exit 1
fi
if [ -z "${flowcell+x}" ]
    then printf "%s\n\nEmply value for -f\n" "$usage" >&2; exit 1
fi
if [ -z "${kit+x}" ]
    then printf "%s\n\nEmply value for -k\n" "$usage" >&2; exit 1
fi
if [ -z "${output+x}" ]
    then printf "%s\n\nEmply value for -o\n" "$usage" >&2; exit 1
fi
echo "Checking parameters: OK"

echo "Creating output directory"
mkdir -p "$output"
echo "Creating output directory: OK"

echo "Writing the qsub script"
uuid=$(uuidgen)
cat <<- EOF > "$output/$uuid.sh"
    #!/usr/bin/env bash
    module load albacore
    read_fast5_basecaller.py -i ${directory} -t ${cpus} -s ${output} \
        -f ${flowcell} -k ${kit} -r -o fastq
EOF
echo "Writing the qsub script: OK"

echo "Queuing albacore"
qsub -N albacore -l h_rt="$run_time":0:0,h_vmem="$mem"G -pe smp "$cpus" -cwd \
    -j y -m sea -M "$email" "$output/$uuid.sh"
echo "Queuing albacore: OK"
exit 0
