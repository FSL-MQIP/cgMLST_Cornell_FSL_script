#!/bin/bash 

set -a # all variables will be exported


### SETTINGS ### 

# Create ID for run; this will just be the current UNIX time. Assumes we will
# never start two runs in the same second.
id=$(date +'%s')
date="$(date)"

# path to genome FASTA is provided as argument
genome=$1

#Location of relevant directories and files. ****YOU MUST CREATE THIS DIRECTORY AHEAD OF TIME AND ADD THE FULL PATH TO IT HERE****
work="/PATH_TO/cgMLST"
dbs=$work/dbs
processed_dir="$work"/processed
#tmp="$work"/tmp
results="$processed_dir"/"$id"_"$(basename $genome)".results
log="$processed_dir"/"$id"_"$(basename $genome)".log
blastdb_log="$processed_dir"/"$id"_"$(basename $genome)".blastdb

# Store this stuff in shared memory to limit I/O
# Check if the directory exists first
#if [ ! -d /dev/shm/lmo-mlst ]
#then
#  mkdir /dev/shm/lmo-mlst
#  ln -s /dev/shm/lmo-mlst "$tmp"
#fi

# Make processed dir if we need to
mkdir -p "$processed_dir"

# Begin log and touch results
printf "Genome:\t"$genome"\nStart:\t$date\n" > "$log"
touch $results

# Number of threads to use for BLAST
num_threads=1

## Copy our files to tmp
#printf "Copying files to memory...\n"
#cp "$genome" "$tmp"/
#cp -r "$dbs" "$tmp"/
#dbs="$tmp"/"$(basename dbs)"
g_base="$(basename $genome)"
#genome="$tmp"/"$(basename $genome)"

### FUNCTIONS ###
# Stuff that's used later but is and hard to read in place

refs() {
  # When adding a new allele, we are going to use a reference of each
  # unique length that occurs for the gene. This function takes the gene
  # name and returns an (arbitrary) example sequence from the database
  # for each length that occurs. These are in descending order of length.
  g_seqs="$(cat $dbs/$1/*.fas)"
  lengths=($(grep '>' <<<"$g_seqs" | cut -f3 -d_ | sort -r | uniq))
  for l in "${lengths[@]}"
  do
    seq="$(grep -A 1 "_$l" <<<"$g_seqs" | head -n2 | tail -n1)"
    printf "$seq\n"
  done
}

match_length_pct() {
  # Takes a line of BLAST output and returns the percent length of the match
  match="$1"
  slength="$(cut -f2 <<<"$match" | cut -f3 -d_)"
  mlength="$(cut -f4 <<<"$match")"
  pct="$( bc -l <<<" "$mlength"/"$slength" " )"
  printf "$pct"
}

match_allele() {
  # Takes a line of BLAST output and returns the allele # of the match
  a="$(cut -f2 <<<"$1" | cut -f2 -d_)"
  printf "$a"
}

match_query() {
  # Takes a line of BLAST output and returns the matched query sequence,
  # reverse complemented if needed
  sstart="$(cut -f9 <<<"$1")"
  send="$(cut -f10 <<<"$1")"
  # Is this contig reverse complemented?
  rev_comp="$(( $sstart > $send ))"

  if [ $rev_comp -eq 1 ]
  then
    q="$(cut -f13 <<<"$1" | tr ACGT TGCA | rev)"
  else
    q="$(cut -f13 <<<"$1")"
  fi

  printf "$q"
}

match_query_plus() {
  # Takes a line of BLAST output and returns the matched query sequnece +18nt
  # on each end. We have to grab it from the genome FASTA file and reverse
  # complement it if needed.

  taxon="$(cut -f1 <<<"$1")"
  qstart="$(cut -f7 <<<"$1")"
  qend="$(cut -f8 <<<"$1")"
  sstart="$(cut -f9 <<<"$1")"
  send="$(cut -f10 <<<"$1")"

  # Is this contig reverse complemented?
  rev_comp="$(( $sstart > $send ))"

  # Grab contig from FASTA
  contig="$(cat "$genome" | \
    sed -e 's/\(^>.*$\)/#\1#/' | \
    tr -d "\r" | \
    tr -d "\n" | \
    sed -e 's/$/#/' | \
    tr "#" "\n" | \
    sed -e '/^$/d' | \
    grep -A 1 "$taxon" | \
    tail -n1)"

  # Try to go back 18 characters. If we go past the beginning of the contig, go
  # forward one contig at a time until we're in the contig. The extra "-1" is
  # because bash string indexing starts at 0, but the BLAST qstart starts at 1.
  start="$(( $qstart - 1 - 18 ))"
  while [ $start -lt 0 ]
  do
    start="$(( $start + 3 ))"
  done
  # the length of the string will be +18, since we're trying to add 18 on the
  # end. Do something similar to the above to make sure we haven't gone off the
  # end of the contig.
  length="$(( $qend - $start + 18 ))"
  while [ "$(( start + length ))" -gt ${#contig} ]
  do
    length="$(( length - 3 ))"
  done

  match_plus="${contig:start:length}"

  # reverse complement if needed
  if [ $rev_comp -eq 1 ]
  then
    match_plus="$(tr ACGT TGCA <<< "$match_plus" | rev)"
  fi

  printf "$match_plus"

}

start_stop_codons() {
  # Reads through first argument; prints if it runs into a START or STOP
  # codon and outputs the position from the beginning (for start) and end
  # (for stop) in the second column
  seq="$1"
  for i in $(seq 0 3 "$(( $(expr length "$seq") - 3 ))")
  do
    nucl="${seq:$i:3}"
    if [ $nucl == "TAA" ]
    then
      printf "STOP\t$(( $(expr length "$seq") - $i - 3 ))\tfrom end\n"
    elif [ $nucl == "TAG" ]
    then
      printf "STOP\t$(( $(expr length "$seq") - $i - 3 ))\tfrom end\n"
    elif [ $nucl == "TGA" ]
    then
      printf "STOP\t$(( $(expr length "$seq") - $i - 3 ))\tfrom end\n"
    elif [ $nucl == "ATG" ]
    then
      printf "START\t$i\tfrom beginning\n"
    elif [ $nucl == "GTG" ]
    then
      printf "START\t$i\tfrom beginning\n"
    elif [ $nucl == "TTG" ]
    then
      printf "START\t$i\tfrom beginning\n"
    elif [ ! -z "${nucl//[-ATCG]}" ]
    then
      printf "NONGATC\t$i\tfrom beginning\n"
    fi
  done
}

ref_align() {
  ### Takes reference sequence as first argument, "new" sequence as
  ### second argument. If alignment shows a frameshift mutation
  ### (non-multiple-of-three set of gap characters in either reference or new
  ### sequence), then it returns "FRAMESHIFT". Otherwise returns the new
  ### sequence after being aligned to reference sequence, trimmed to its
  ### length, and ungapped.
  ref_seq=">ref\n$1"
  new_seq=">new\n$2"
  # Align reference and new allele
  aligned_seqs="$(printf "$ref_seq\n$new_seq" | muscle -quiet | sed -e\
    's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr\
    "#" "\n" | sed -e '/^$/d')"
  # Remove leading and trailing gaps from reference and save
  ref_trimmed="$(printf "$aligned_seqs" | sed '2q;d' | sed 's/^-\+//' | sed 's/-\+$//')"
  # Grab just leading gaps from ref sequence
  ref_leading_gaps="$(printf "$aligned_seqs" | sed '2q;d' | sed 's/^\(-*\).\+$/\1/')"
  # Grab untrimmed new sequence
  new_untrimmed="$(printf "$aligned_seqs" | sed '4q;d')"
  # Grab only portion of new sequence which does not extend beyond the ref sequence
  new_trimmed="${new_untrimmed:${#ref_leading_gaps}:${#ref_trimmed}}"

  frameshift="$(printf "$ref_trimmed\n$new_trimmed" | sed 's/---//g' | grep '-')"

  if [ ! -z "$frameshift" ]
  then
    printf "FRAMESHIFT"
  elif [ "$3" = "trimmed" ]
  then
    printf "$(printf "$new_trimmed" | sed 's/-//g')"
  else
    printf "$(printf "$new_untrimmed" | sed 's/-//g')"
  fi
}

blaster() {
  gene="$(basename $1)"
  gene_dir="$work"/dbs/"$gene" # This is in the ORIGINAL work directory, not the tmp directory
  # gene_db="$1"/"$gene" # I don't think I need this line
  gene_db="$gene_dir"/"$gene"
  # Get match with highest bitscore

  # Create debug directory            
  debug="$work"/debug/"$id"/"$gene"   
  mkdir -p "$debug"                   

  # First, search for longest 100% match
  # This CRITICALLY REQUIRES the taxon names in the allele database to have
  # _exactly three_ fields, separated by underscores, and that the length of
  # the allele be the third field.
  match="$(blastn -db "$gene_db" \
    -query "$genome" \
    -word_size 10 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" \
    -perc_identity 100 \
    -num_threads 1 | \
    awk '{ split($2,a,"_"); if ( a[3]==$4 ) print $0}' | \
    head -n 1)"
  # Check % of match length
  pct_length="$(match_length_pct "$match")"

  printf "### TOP 100 PCT MATCH:\n$match \n### PCT LENGTH: \n$pct_length\n" > "$debug"/"1. 100pct.txt" # dbg dk657

  # OLD: If the longest 100% id match isn't 100% length, grab the match with the
  # OLD: highest bitscore

  # Revised: If there is no 100% id match, grab the match with the highest
  # bitscore
  #if [ "$(bc -l <<<"$pct_length == 1")" -eq 0 ]
  if [ -z "$match" ]
  then
    match="$(blastn -db "$gene_db" \
      -query "$genome" \
      -word_size 10 \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" \
      -perc_identity 0.7 \
      -num_threads 1 \
      -max_target_seqs 1 | \
      sort -n -k12 -r | \
      head -n1)"
    pct_length="$(match_length_pct "$match")"
  fi

  pct_id="$(cut -f3 <<<"$match")"

  printf "### TOP MATCH (EITHER 100 PCT ID/LENGTH OR HIGHEST BITSCORE):\n$match \n### PCT LENGTH: \n$pct_length \n### PCT ID: \n$pct_id\n" > "$debug"/"2. top match.txt" # dbg dk657

  if [ "$(bc -l <<<"$pct_length == 1")" -eq 1 -a "$(bc -l <<<"$pct_id == 100")" -eq 1 ]
  then
    printf "Exact match found\n" >> "$debug"/process.txt  # dbg dk657
    # Exact match
    allele="$(match_allele "$match")"
    printf "$gene:\tExact match of allele $allele\n"
    printf "$g_base\t$gene\t$allele\n" >> "$results"
  elif [ "$(bc -l <<<"$pct_length > 0.7")" -eq 1 -a "$(bc -l <<<"$pct_id > 70")" -eq 1 ]
  then
    printf "New allele process started.\n" >> "$debug"/process.txt  # dbg dk657
    # New allele
    pass=0 # indicator to skip extended check if we pass using the 

    # Create array with arbitrary references from this gene's db;
    # one for each unique length
    readarray references <<< "$(refs $gene)"
    printf '%s\n' "${references[@]}" >> "$debug"/new_allele_references.txt    # dbg dk657

    # First we'll grab the query sequence from the BLAST match
    query="$(match_query "$match")"

    # Now we'll check that against one reference of each length
    for r in "${references[@]}"
    do
      r="$(printf "$r" | tr -d "\n" | tr -d "\r")" # delete newline character >:(
      # Frameshift/get seq
      seq_aligned="$(ref_align "$r" "$query" trimmed)"
      seq_aligned="$(printf "$seq_aligned" | tr -d "\n" | tr -d "\r")"
      start_stop="$(start_stop_codons "$seq_aligned")"
      start="$(grep START <<<"$start_stop")"
      stop="$(grep STOP <<<"$start_stop")"
      nongatc="$(grep NONGATC <<<"$start_stop")"

      printf ">ref_"${#r}"\n$r\n>$gene\n$query\n" > "$debug"/ref_"${#r}".fa   # dbg dk657
      printf "$start_stop" > "$debug"/ref_"${#r}"_startstop.txt               # dbg dk657
      printf "$seq_aligned" > "$debug"/ref_"${#r}"_seq_aligned.txt            # dbg dk657

      if [ "$seq_aligned" = "FRAMESHIFT" ]
      then
        # FRAMESHIFT
        message="With "${#r}"nt reference, sequence for "$gene" cannot be added because there is a frameshift mutation.\n"
        printf "$message" >> "$log"
      elif [ -z "$start" ]
      then
        # NO START CODON AT ALL
        message="With "${#r}"nt reference, sequence for "$gene" cannot be added because there is no start codon.\n"
        printf "$message" >> "$log"
      elif [ -z "$stop" ]
      then
        # NO STOP CODON AT ALL
        message="With "${#r}"nt reference, sequence for "$gene" cannot be added because there is no stop codon.\n"
        printf "$message" >> "$log"
      elif [ "$(sort -k2 -n <<<"$start" | head -n1 | cut -f 2)" -gt 15 ]
      then
        # NO START CODON IN FIRST SIX
        message="With "${#r}"nt reference, sequence for "$gene" cannot be added because there is no start codon present in the first six codons.\n"
        printf "$message" >> "$log"
      elif [ "$(sort -k2 -n -r <<<"$stop" | head -n1 | cut -f 2)" -gt 15 ]
      then
        # PREMATURE STOP CODON
        message="With "${#r}"nt reference, sequence for "$gene" cannot be added because there is a stop codon more than six codons from the end.\n"
        printf "$message" >> "$log"
      elif [ ! -z "$nongatc" ]
      then
        # NON-GATC CHARACTERS
        message="With "${#r}"nt reference, sequence for "$gene" cannot be added because there are non-GATC characters in the sequence.\n"
        printf "$message" >> "$log"
      else
        # SUCCESS
        # We need to 1. Mark that the allele passed so
        # that we don't grab the exteded sequence, 2.
        # add the sequence to the database, and 3.
        # break the loop so we don't check against any
        # of the other references.
        pass=1

        # Last allele in db
        last_allele="$(grep '>' <"$gene_db".fas | sed \
          's/_/\t/g' | sort -k2 -n -r | head -n1 | cut -f \
          2)" 
        # New allele #
        new_allele="$(( $last_allele + 1 ))" 
        # take START closest to beginning
        start_pos="$(sort -k2 -n <<<"$start" | head -n1 \
          | cut -f 2)" 
        # take STOP farthest from end
        end_pos="$(( $(expr length "$seq_aligned") - \
          $(sort -k2 -n -r <<<"$stop" | head -n1 | cut -f \
          2) - $start_pos ))" 
        # Sequence for new allele
        new_allele_seq="$( sed 's/-//g' \
          <<<"${seq_aligned: $start_pos : $end_pos }" )"
        new_allele_length="$(expr length \
          $new_allele_seq)"

        # Write new allele to a FASTA file in the gene's dir
        printf ">$gene"_"$new_allele"_"$new_allele_length\n""$new_allele_seq\n" > "$gene_dir"/chunk_"$id".fas

        # Re-build database
        cat "$gene_dir"/original_"$gene".fas "$gene_dir"/chunk_* > "$gene_dir"/"$gene".fas
        makeblastdb -in "$gene_dir"/"$gene".fas -dbtype nucl -title "$gene" -out "$gene_dir"/"$gene" >> "$blastdb_log"

        message="With "${#r}"nt reference, sequence for "$gene" has been added as "$new_allele_length"nt new allele "$new_allele".\n\n"
        printf "$message" >> "$log"
        printf "$message"

        printf "$g_base\t$gene\t$new_allele\n" >> "$results"
        break;
      fi
    done

    # Now, if the new allele didn't pass against any of the
    # references, we'll grab the longer version (+18nt on either
    # end) and check again.

    # Read in extended match query
    query_extended="$(match_query_plus "$match")"

    # Again, we'll check that against one reference of each length
    # This code is highly repetitive and needs to be put in a function.
    for r in "${references[@]}"
    do
      if [ "$pass" -eq 1 ]
      then
        # If we already passed without the extended
        # sequence, don't do any of this
        break;
      fi
      r="$(printf "$r" | tr -d "\n" | tr -d "\r")" # delete newline character >:(
      # Frameshift/get seq
      seq_aligned="$(ref_align "$r" "$query_extended" untrimmed)"
      seq_aligned="$(printf "$seq_aligned" | tr -d "\n" | tr -d "\r")"
      start_stop="$(start_stop_codons "$seq_aligned")"
      start="$(grep START <<<"$start_stop")"
      stop="$(grep STOP <<<"$start_stop")"
      nongatc="$(grep NONGATC <<<"$start_stop")"
      if [ "$seq_aligned" = "FRAMESHIFT" ]
      then
        # FRAMESHIFT
        message="With "${#r}"nt reference, extended sequence for "$gene" cannot be added because there is a frameshift mutation.\n"
        printf "$message" >> "$log"
        printf "$message"
      elif [ -z "$start" ]
      then
        # NO START CODON AT ALL
        message="With "${#r}"nt reference, extended sequence for "$gene" cannot be added because there is no start codon.\n"
        printf "$message" >> "$log"
        printf "$message"
      elif [ -z "$stop" ]
      then
        # NO STOP CODON AT ALL
        message="With "${#r}"nt reference, extended sequence for "$gene" cannot be added because there is no stop codon.\n"
        printf "$message" >> "$log"
        printf "$message"
      elif [ "$(sort -k2 -n <<<"$start" | head -n1 | cut -f 2)" -gt 15 ]
      then
        # NO START CODON IN FIRST SIX
        message="With "${#r}"nt reference, extended sequence for "$gene" cannot be added because there is no start codon present in the first six codons.\n"
        printf "$message" >> "$log"
        printf "$message"
      elif [ "$(sort -k2 -n -r <<<"$stop" | head -n1 | cut -f 2)" -gt 15 ]
      then
        # PREMATURE STOP CODON
        message="With "${#r}"nt reference, extended sequence for "$gene" cannot be added because there is a stop codon more than six codons from the end.\n"
        printf "$message" >> "$log"
        printf "$message"
      elif [ ! -z "$nongatc" ]
      then
        # NON-GATC CHARACTERS
        message="With "${#r}"nt reference, extended sequence for "$gene" cannot be added because there are non-GATC characters in the sequence.\n"
        printf "$message" >> "$log"
        printf "$message"
      else
        # SUCCESS
        # We need to 1. Mark that the allele passed so that we don't grab the exteded sequence, 2.  add the sequence to the database, and 3.  break the loop so we don't check against any of the other references.
        pass=1

        # Last allele in db
        last_allele="$(grep '>' <"$gene_db".fas | sed \
          's/_/\t/g' | sort -k2 -n -r | head -n1 | cut -f \
          2)" 
        # New allele #
        new_allele="$(( $last_allele + 1 ))" 
        # take START closest to beginning
        start_pos="$(sort -k2 -n <<<"$start" | head -n1 \
          | cut -f 2)" 
        # take STOP farthest from end
        end_pos="$(( $(expr length "$seq_aligned") - \
          $(sort -k2 -n -r <<<"$stop" | head -n1 | cut -f \
          2) - $start_pos ))" 
        # Sequence for new allele
        new_allele_seq="$( sed 's/-//g' \
          <<<"${seq_aligned: $start_pos : $end_pos }" )"
        new_allele_length="$(expr length \
          $new_allele_seq)"

        # Write new allele to a FASTA file in the gene's dir
        printf ">$gene"_"$new_allele"_"$new_allele_length\n""$new_allele_seq\n" > "$gene_dir"/chunk_"$id".fas

        # Re-build database
        cat "$gene_dir"/original_"$gene".fas "$gene_dir"/chunk_* > "$gene_dir"/"$gene".fas
        makeblastdb -in "$gene_dir"/"$gene".fas -dbtype nucl -title "$gene" -out "$gene_dir"/"$gene" >> "$blastdb_log"

        message="With "${#r}"nt reference, extended sequence for "$gene" has been added as "$new_allele_length"nt new allele "$new_allele".\n\n"
        printf "$message" >> "$log"
        printf "$message"
        printf "$g_base\t$gene\t$new_allele\n" >> "$results"
        break;
      fi
    done

    if [ "$pass" -eq 0 ]
    then
      message="Sequence for "$gene" found but did not pass criteria to be added as new allele.\n\n"
      printf "$message" >> "$log"
      printf "$message"
      printf "$g_base\t$gene\tNA\n" >> "$results"
    fi
  else
    # No match
    printf "$gene:\tNot found\n"
    printf "$g_base\t$gene\tNA\n" >> "$results"
    printf "No match found for $gene\n\n" >> "$log"
  fi 
}


set +a

### RUN IT ##

printf "BLASTing genome...\n"
genes=($dbs/*)
parallel -P $num_threads blaster ::: "${genes[@]}"

# Final message for log:
new_alleles="$(grep "has been added" "$log" | sed 's/^.\+\(lmo[0-9]\{4\}\).\+$/\1/' | sort)"
not_added="$(grep "did not pass" "$log" | sed 's/^.\+\(lmo[0-9]\{4\}\).\+$/\1/' | sort)"
missing="$(grep "No match" "$log" | sed 's/^.\+\(lmo[0-9]\{4\}\)$/\1/' | sort)"

printf "New alleles were found and added to the database for the following "$(printf "$new_alleles" | wc -l)" genes.\n" >> "$log"
printf "$new_alleles\n" >> "$log"

printf "Matches were found for the following "$(printf "$not_added" | wc -l)" genes, but the sequences did not meet criteria for addition to the database:\n" >> "$log"
printf "$not_added\n" >> "$log"

printf "Matches were not found for the following "$(printf "$missing" | wc -l)" genes:\n" >> "$log"
printf "$missing\n" >> "$log"
