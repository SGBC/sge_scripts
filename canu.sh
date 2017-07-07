#!/usr/bin/env bash

# assembly for long reads (using canu)

usage="$(basename "$0") [-h] [-e email] [-c cpus] [-m ram] [-t time]
    [-i reads.fastq] [-t read type] [-o output] [-g genome size] --

        assembly for long reads (using canu)

        qsub options:
    		-h show this useful help
            -e <string> your email address
            -c <int> number of cpus to use
            -m <int> number of ram per cpu to use
            -t <time> upper time limit for the job, in hours

        albacore options:
    		-i <fastq> input fastq file
            -r <string> read type. Can be pacbio-raw, pacbio-corrected,
                nanopore-raw or nanopore-corrected
            -g <int> approx size of genome
            -o <path> path to save output"

while getopts "he:c:m:t:i:r:g:o:" option
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
            total_mem=$(echo $cpus*$mem | bc -l )G
                ;;
            t) run_time="$OPTARG"
                ;;
            i) reads=$(readlink -f "$OPTARG")
            prefix_reads=$(basename ${reads%.*})
                ;;
            r) read_type="pacbio-raw pacbio-corrected \
            nanopore-raw nanopore-corrected"
            [[ $read_type =~ $OPTARG ]] && type="$OPTARG" || \
            { printf "%s\n\Read type not recognised\n" "$usage">&2; exit 1; }
                ;;
            g) genome_size="$OPTARG"
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
if [ -z "${reads+x}" ]
    then printf "%s\n\nEmply value for -i\n" "$usage" >&2; exit 1
fi
if [ -z "${type+x}" ]
    then printf "%s\n\nEmply value for -r\n" "$usage" >&2; exit 1
fi
if [ -z "${genome_size+x}" ]
    then printf "%s\n\nEmply value for -g\n" "$usage" >&2; exit 1
fi
if [ -z "${output+x}" ]
    then printf "%s\n\nEmply value for -o\n" "$usage" >&2; exit 1
fi
echo "Checking parameters: OK"

echo "Creating output directory"
mkdir -p "$output" || exit 1
echo "Creating output directory: OK"

echo "Writing the qsub script"
uuid=$(uuidgen)
cat <<- EOF > "$output/$uuid.sh"
    #!/usr/bin/env bash
    module load canu
    canu -p ${prefix_reads} -d ${output} \
        genomeSize=${genome_size} -${read_type} ${reads} \
        maxThreads=${cpus} maxMemory=${total_mem}
EOF
echo "Writing the qsub script: OK"

echo "Queuing canu"
qsub -N canu -l h_rt="$run_time":0:0,h_vmem="$mem"G -pe smp "$cpus" -cwd \
    -j y -m sea -M "$email" "$output/$uuid.sh"
echo "Queuing canu: OK"
exit 0
