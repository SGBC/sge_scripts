#!/usr/bin/env bash

# bowtie2 for paired-end reads

usage="$(basename "$0") [-h] [-e email] [-c cpus] [-m ram] [-t time]
    [-r reference] [-1 reads_1.fq] [-2 reads_2.fq] [-o output] --

        bowtie2 for paired-end reads

        qsub options:
    		-h show this useful help
            -e <string> your email address
            -c <int> number of cpus to use
            -m <int> number of ram per cpu to use
            -t <time> upper time limit for the job, in hours

        bowtie options:
    		-r <fasta> reference in fasta format
            -1 <fastq> forward reads
            -2 <fastq> reverse reads (optional, only use -1 for single end)
    		-o <path> output directory"

while getopts "he:c:m:t:r:1:2o:" option
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
            # tot_mem=$(echo $cpus*$mem | bc -l )G
                ;;
            t) run_time="$OPTARG"
                ;;
            r) reference=$(readlink -f "$OPTARG")
    		prefix_ref=$(basename ${reference%.*})
                ;;
            1) r1=$(readlink -f "$OPTARG")
            prefix_reads=$(basename ${r1%.*})
                ;;
            2) r2=$(readlink -f "$OPTARG")
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
if [ -z "${r1+x}" ]
    then printf "%s\n\nEmply value for -1\n" "$usage" >&2; exit 1
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
if [ -z "${r2+x}" ]
    then
        cat <<- EOF > "$output/$uuid.sh"
        #!/usr/bin/env bash
        module load bowtie
        module load samtools
        if [ ! -f "$output/$prefix_ref".1.bt2 ]
            then
                bowtie2-build "$reference" "$output/$prefix_ref"
        fi
        bowtie2 -p "$cpus" -x "$output/$prefix_ref" -U "$r1" | \
        samtools view -@ "$cpus" -b - | samtools sort -@ "$cpus" -m "$mem"G \
            -o "$output/$prefix_reads.bam"
EOF
    else
        cat <<- EOF > "$output/$uuid.sh"
        #!/usr/bin/env bash
        module load bowtie
        module load samtools
        if [ ! -f "$output/$prefix_ref".1.bt2 ]
            then
                bowtie2-build "$reference" "$output/$prefix_ref"
        fi
        bowtie2 -p "$cpus" -x "$output/$prefix_ref" -1 "$r1" -2 "$r2" | \
        samtools view -@ "$cpus" -b - | samtools sort -@ "$cpus" -m "$mem"G \
            -o "$output/$prefix_reads.bam"
EOF

fi

echo "Writing the qsub script: OK"

echo "Queuing Bowtie"
qsub -N Bowtie -l h_rt="$run_time":0:0,h_vmem="$mem"G -pe smp "$cpus" -cwd \
    -j y -m sea -M "$email" "$output/$uuid.sh"
echo "Queuing Bowtie: OK"
exit 0
