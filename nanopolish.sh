#!/usr/bin/env bash

# consensus step for OLC ONT assembly

usage="$(basename "$0") [-h] [-e email] [-c cpus] [-m ram] [-t time]
    [-i reads.fastq] [-t read type] [-o output] [-g genome size] --

        consensus step for OLC ONT assembly (using nanopolish)

        qsub options:
    		-h show this useful help
            -e <string> your email address
            -c <int> number of cpus to use
            -m <int> number of ram per cpu to use
            -t <time> upper time limit for the job, in hours

        nanopolish options:
    		-i <fastq> input fastq|a file
            -a <fasta> assembly in fasta format
            -o <path> path to save output"

while getopts "he:c:m:t:i:a:o:" option
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
            i) reads=$(readlink -f "$OPTARG")
                ;;
            a) assembly=$(readlink -f "$OPTARG")
            prefix_assembly=$(basename ${assembly%.*})
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
if [ -z "${assembly+x}" ]
    then printf "%s\n\nEmply value for -a\n" "$usage" >&2; exit 1
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
    module load nanopolish
    module load minimap/2.1
    module load samtools
    module load parallel
    minimap -ax map-ont -t "${cpus}" "${assembly}" "${reads}" | \
    samtools view -@ "${cpus}" -b - | samtools sort -@ "${cpus}" -m "${mem}G" \
        -o "${output}/${prefix_assembly}.bam"
    samtools index "${output}/${prefix_assembly}.bam"
    nanopolish_makerange.py "${assembly}" | \
    parallel --results "${output}" -P "${cpus}" \
    nanopolish variants --consensus "${output}/"polished.{1}.fa -w {1} -r "${reads}" \
        -b "${output}/${prefix_assembly}.bam" -g "${assembly}" -t "${cpus}" \
        --min-candidate-frequency 0.1
    nanopolish_merge.py "${output}/polished.*.fa" > \
        "${output}/${prefix_assembly}_polished.fasta"
EOF
echo "Writing the qsub script: OK"

echo "Queuing nanopolish"
qsub -N nanopolish -l h_rt="$run_time":0:0,h_vmem="$mem"G -pe smp "$cpus" -cwd \
    -j y -m sea -M "$email" "$output/$uuid.sh"
echo "Queuing nanopolish: OK"
exit 0
