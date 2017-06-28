#!/usr/bin/env bash

# long read mapping (with graphmap)

usage="$(basename "$0") [-h] [-e email] [-c cpus] [-m ram] [-t time]
    [-r reference] [-1 reads_1.fq] [-2 reads_2.fq] [-o output] --

        long read mapping (with graphmap)

        qsub options:
    		-h show this useful help
            -e <string> your email address
            -c <int> number of cpus to use
            -m <int> number of ram per cpu to use
            -t <time> upper time limit for the job, in hours

        graphmap options:
    		-r <fasta> reference genome
            -d <fastq> path to the reads file.
            -o <path> path to save output"

while getopts "he:c:m:t:r:d:o:" option
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
            r) reference=$(readlink -f "$OPTARG")
                ;;
            d) reads=$(readlink -f "$OPTARG")
            prefix_reads=$(basename ${reads%.*})
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
if [ -z "${reference+x}" ]
    then printf "%s\n\nEmply value for -r\n" "$usage" >&2; exit 1
fi
if [ -z "${prefix_reads+x}" ]
    then printf "%s\n\nEmply value for -d\n" "$usage" >&2; exit 1
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
    module load graphmap
    graphmap align -r ${reference} -t ${cpus} -d ${reads} \
        -o ${output}/${prefix_reads}.sam
EOF
echo "Writing the qsub script: OK"

echo "Queuing graphmap"
qsub -N graphmap -l h_rt="$run_time":0:0,h_vmem="$mem"G -pe smp "$cpus" -cwd \
    -j y -m sea -M "$email" "$output/$uuid.sh"
echo "Queuing graphmap: OK"
exit 0
