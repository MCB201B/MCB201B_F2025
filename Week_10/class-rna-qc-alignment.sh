#!/bin/bash

# This is a bash script that will run through QC and alignment.

# Confirm location
output_parent=~/MCB201B_F2025/Week_10
echo -e "Changing to this directory: ${output_parent} \n"
cd "${output_parent}"

# Make directories as needed to hold our outputs
echo -e "Making directories to keep things tidy...\n"
mkdir -p "${output_parent}/seq-view"
mkdir -p "${output_parent}/fastqc-results"

mkdir -p "${output_parent}/alignments-temp"
mkdir -p "${output_parent}/alignment-logs"

alignment_output="${output_parent}/alignments-temp"

mkdir -p "${alignment_output}/sams"
mkdir -p "${alignment_output}/bams"

echo -e "Setting our file path to our RNA-seq data...\n"
file_path=~/shared/2025-fall/courses/1547808/rna-seq/truncated/g1
echo -e "Confirming our file path as: \n   ${file_path}"
cd "${file_path}"

# Pull our file names into arrays
mate_1=(*_r1.fastq.gz)
mate_2=(*_r2.fastq.gz)
full_set=(*.fastq.gz)

# Check to see if we successfully pulled in files
# If we did, then variable should not be an empty string
# So -z which assesses whether or not something is an empty string,
# will not be satisfied, and script will proceed.
if [[ -z "${mate_1[@]}" ]]; then
    echo "Oh my god Help!"
    
    # Force our script to exit with value 1, indicating failure.
    exit 1
else
    echo -e "Pulled in the following files to work on:"
    echo -e "  ${mate_1[@]}"
    echo -e "  ${mate_2[@]}"
fi

# Confirm the index
index_dir=~/shared/2025-fall/courses/1547808/rna-seq/rna-index/hg19
echo "Confirming index directory and base names as:"
echo "  ${index_dir}"

# This is an additional add on from end of Friday
# This conditional statement is set up to skip alignment if we have bams
# You can see that it's currently set for the truncated dataset, so you'll need to change.
if [[ -f "${alignment_output}/bams/1M_g1_ctrl-name.bam" \
    && -f "${alignment_output}/bams/1M_g1_tazko-name.bam" ]]; then
    
    # Our htseq-count later on uses the variable base_name
    # So we can set it up here, so we don't run into issues later
    # Since we would otherwise skip the base_name assignment part of our script
    base_name=1M_g1
    
else
    # Operate on individual files first
    for file_name in "${full_set[@]}"; do
        echo "Pulling out ten reads from..."
        base_file_name=${file_name%.fastq.gz}
        echo "  ${base_file_name}"
        zcat "${file_name}" \
            | head -40 \
            > "${output_parent}/seq-view/${base_file_name}-10-reads.txt"
    
        echo -e "Running a FastQC for ${base_file_name}...\n"
        fastqc \
            -o "${output_parent}/fastqc-results" \
            "${file_name}"
    done
    
    # Initialize our for loop for alignment and playing with samtools
    for ((i=0; i<"${#mate_1[@]}"; i++)); do
        # Obtain our new basename
        base_name="${mate_1[i]%_r1.fastq.gz}"
        echo -e "Running alignment on ${base_name}...\n"
    
        hisat2 \
            -x "${index_dir}" \
            -1 "${mate_1[i]}" \
            -2 "${mate_2[i]}" \
            --rna-strandness RF \
            -p 8 \
            --mm \
            -S "${alignment_output}/sams/${base_name}.sam" \
            --summary-file "${output_parent}/alignment-logs/${base_name}-alignment.log" \
            --new-summary

        # Here, we set up a checkpoint to confirm if alignment was successfully completed
        # We do so by checking for the existence of our alignment log file using -f
        # -f simply checks if it is a file, and can be used to check if a file exists or not
        # hisat2 will not create any log if alignment did not complete
        if [[ -f "${output_parent}/alignment-logs/${base_name}-alignment.log" ]]; then
            echo -e "Alignment for ${base_name} has been completed.\n"
        else
            echo -e "Alignment for ${base_name} has an issue! ABORT! \n"
            exit 1
        fi
    
        # Convert SAM to a BAM
        samtools view \
            -b \
            -@ 8 \
            -o "${alignment_output}/bams/${base_name}.bam" \
            "${alignment_output}/sams/${base_name}.sam"
    
        # Pull lines out of SAM to explore later
        head \
            -150 \
            "${alignment_output}/sams/${base_name}.sam" \
            > "${output_parent}/alignment-logs/${base_name}-150-rows-SAM.txt"
    
        # Sort our BAM by position
        samtools sort \
            -@ 4 \
            -o "${alignment_output}/bams/${base_name}-position.bam" \
            "${alignment_output}/bams/${base_name}.bam"

        # Set up checkpoint to confirm BAM was sorted, which means no issues
        # So then we can go ahead and delete our SAM and proceed with rest of script
        if [[ -f "${alignment_output}/bams/${base_name}-position.bam" ]]; then
            rm -f "${alignment_output}/sams/${base_name}.sam"
            echo -e "Successfully sorted BAM so now deleted SAM for ${base_name} \n"
        else
            echo -e "Issue with BAM file somewhere! ABORT! \n"
            exit 1
        fi
    
        # Index our position-sorted BAM
        samtools index \
            "${alignment_output}/bams/${base_name}-position.bam"
    
        # Sort our BAM by name for later read counting
        samtools sort \
            -n \
            -@ 4 \
            -o "${alignment_output}/bams/${base_name}-name.bam" \
            "${alignment_output}/bams/${base_name}.bam"

        # Set up another checkpoint
        # This time THREE conditions must be satisfied
        # Indicated by the double ampersand indicating AND
        # We are checking if we have our sorted BAMs and indexed the position-sorted one
        # If so, we can delete our unsorted BAM and proceed with rest of script
        if [[ -f "${alignment_output}/bams/${base_name}-position.bam" \
            && -f "${alignment_output}/bams/${base_name}-position.bam.bai" \
            && -f "${alignment_output}/bams/${base_name}-name.bam" ]]; then
                echo -e "Everything looks good! \n"
                rm -f "${alignment_output}/bams/${base_name}.bam"
        else
            echo -e "Something went horribly wrong! ABORT! \n"
            exit 1
        fi
    
        echo -e "Alignment is completed along with BAM sorting! \n"
    done
fi

# Verify our HTSeq installation
# We can set up a checkpoint to verify installation
# Recall that if a command errors out it exits with a value that is NOT zero
# 0 indicates success
# Anything 1 or above indicates failure
# We can take advantage of this to check to see if we have HTSeq installed
# by calling up the exit value of the last command with "${?}"
# basically the special character variable that holds exit value of the most recent command that was run
# So we can set up a conditional to see if it equals 0 or not when we run htseq-count to check what version we have

# So first run our command, which will produce some exit value
htseq-count --version

# Then we can perform our check using the exit value
if [[ "${?}" == 0 ]]; then

    # If it is 0, then it successfully checked the version, and we know it's installed
    echo -e "HTSeq is installed. Proceeding with counting...\n"
else

    # If it can't, then that means it's not installed, so we proceed with installation.
    echo -e "HTSeq is not yet installed. Proceeding with installation...\n"
    pip install HTSeq

    # Then we can do one more check.
    htseq-count --version
    if [[ "${?}" == 0 ]]; then
        echo -e "HTSeq is installed. Proceeding with counting...\n"
    else
        echo -e "Big UH OH. Can't install HTSeq... ABORT! \n"
        exit 1
    fi
fi

# Now we confirmed HTSeq is installed, and we can proceed with read counting.

# Make additional needed directories
mkdir -p "${output_parent}/counts"

# Pull in GTF file
gtf_file=~/shared/2025-fall/courses/1547808/rna-seq/rna-feature/hg19-refseq.gtf


# This is to output some rows from GTF file to explore in class.
cat "${gtf_file}" \
    | head -100 \
    > "${output_parent}/sample-gtf-to-explore.txt"

# Pull in our name-sorted BAMs
cd "${alignment_output}/bams"
sorted_bams=(*-name.bam)

# Run htseq-count on our two sample groups
htseq-count \
    -t exon \
    -i gene_id \
    -r name \
    -s reverse \
    -f bam \
    "${sorted_bams[@]}" \
    "${gtf_file}" \
    > "${output_parent}/counts/${base_name%_tazko}_counts.txt"

echo -e "Finally finished...\n"
