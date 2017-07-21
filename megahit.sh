#!/usr/bin/env bash

# (meta)genome assembly using megahit

usage="$(basename "$0") [-h] [-e email] [-c cpus] [-m ram] [-t time] --

    (meta)genome assembly using megahit

        qsub options:
            -h show this useful help
            -e <string> your email address
            -c <int> number of cpus to use
            -m <int> number of ram per cpu to use
            -t <time> upper time limit for the job, in hours

        megahit options:
            -1 <fastq> reads file in fastq format. Can be a comma-separated
                list of forward reads if several libraries
            -2 <fastq> reads file in fastq format. Same as for -1
            -o <path> output directory"

while getopts "he:c:m:t:p1:2:o:" option
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
            tot_mem=$(echo $cpus*$mem*1000000000 | bc -l )
                ;;
            t) run_time="$OPTARG"
                ;;
            1) r1="$OPTARG"
                ;;
            2) r2="$OPTARG"
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
if [ -z "${tot_mem+x}" ]
    then printf "%s\n\Bad value for -m\n" "$usage" >&2; exit 1
fi
if [ -z "${run_time+x}" ]
    then printf "%s\n\nEmply value for -t\n" "$usage" >&2; exit 1
fi
if [ -z "${r1+x}" ]
    then printf "%s\n\nEmply value for -1\n" "$usage" >&2; exit 1
fi
if [ -z "${r2+x}" ]
    then printf "%s\n\nEmply value for -2\n" "$usage" >&2; exit 1
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
    module load megahit
    megahit -t "${cpus}" -m "${tot_mem}" -1 "${r1}" -2 "${r2}" \
        -f -o "${output}"
EOF
echo "Writing the qsub script: OK"

echo "Queuing Megahit"
qsub -N Megahit -l h_rt="$run_time":0:0,h_vmem="$mem"G -pe smp "$cpus" -cwd \
    -j y -m sea -M "$email" "$output/$uuid.sh"
echo "Queuing Megahit: OK"
exit 0
