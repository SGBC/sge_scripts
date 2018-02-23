#!/usr/bin/env bash

# basecalling for nanopore reads (with albacore)

usage="$(basename "$0") [-h] [-e email] [-c cpus] [-m ram] [-t time]
    [-l fastq files] [-o output]  --

        Unicycler for nanopore reads

        qsub options:
    		-h show this useful help
            -e <string> your email address
            -c <int> number of cpus to use
            -m <int> number of ram per cpu to use
            -t <time> upper time limit for the job, in hours

        Unicycler options:
    		-l <string> fastq file with trimmed reads
        -o <path> path to save output"

while getopts "he:c:m:t:l:o:" option
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
            l) long_reads=$(readlink -f "$OPTARG")
                ;;
            o) output=$(readlink -f "$OPTARG")
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
if [ -z "${long_reads+x}" ]
    then printf "%s\n\nEmply value for -i\n" "$usage" >&2; exit 1
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
    module load unicycler
    unicycler -l ${long_reads} -t ${cpus} -o ${output} \
    --verbosity 2 --keep 0
EOF
echo "Writing the qsub script: OK"

echo "Queuing Unicycler"
qsub -N Unicycler -l h_rt="$run_time":0:0,h_vmem="$mem"G -pe smp "$cpus" -cwd \
    -j y -m sea -M "$email" "$output/$uuid.sh"
echo "Queuing Unicycler: OK"
exit 0
