#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH -N 4
#SBATCH -n 10
#SBATCH --mem-per-cpu 8g
#SBATCH -J 20250412_21-28-hmmer-dmr
#SBATCH --output jobreports/20250412_21-28-hmmer-dmr-out.txt
#SBATCH --error jobreports/20250412_21-28-hmmer-dmr-out.txt
# I am requesting 4 nodes containing 10 CPUs, with 8 GB memory per CPU. Total: 80 GB

source ~/miniconda3/etc/profile.d/conda.sh
# Install all tools with conda
## Setup channels for conda installation
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

## Create a conda environment
conda create -yqn nfkappab_evo-tools -c conda-forge -c bioconda \
    blast bedtools emboss transdecoder hmmer clustalo

date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
date_time=20250415_18_50_06
nthreads=10
interpro_output_date_time=20250409
full_fasta_date_time=20250409
# transdecoder_date_time=${date_var}_${time_var}
transdecoder_date_time=20250328_13_43_35
hmmer_date_time=${date_var}_${time_var}
hmmer_date_time=20250415_18_50_06
hmmer_date=20250410

project_home_dir=~/projects/metagenome
transdecoder_output_dir=output/transdecoder_output
interpro_fasta_path=output/fasta/"${interpro_output_date_time}"/interpro
full_fasta_path=output/fasta/"${full_fasta_date_time}"
alignments_dir=output/alignments/"${hmmer_date}"
hmmprofiles_dir=output/hmmprofiles/"${hmmer_date}"
hmmstdout_dir=output/hmmstdout/"${hmmer_date}"
hmmtables_dir=output/hmmtables/"${hmmer_date}"
hmmalignments_dir=output/hmmalignments/"${hmmer_date}"
hmmfasta_dir=output/hmmfasta/"${hmmer_date}"
hmmscan_out_dir=output/hmmscan_out/"${hmmer_date}"
ssi_index_dir=output/ssindex/"${date_var}"
interpro_hmms_dir=data/interpro_hmms

mkdir -p data/ncbi_datasets 
mkdir -p jobreports
mkdir -p output/fasta
mkdir -p "${transdecoder_output_dir}" 
mkdir -p "${alignments_dir}"
mkdir -p "${hmmprofiles_dir}"
mkdir -p "${hmmstdout_dir}"
mkdir -p "${hmmtables_dir}"
mkdir -p "${hmmalignments_dir}"
mkdir -p "${hmmscan_out_dir}"
mkdir -p "${ssi_index_dir}"
mkdir -p "${interpro_hmms_dir}"

# Download a pre-built HMM from InterPro Pfam
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz \
 -O "${interpro_hmms_dir}"/Pfam-A.hmm.gz
gunzip "${interpro_hmms_dir}"/Pfam-A.hmm.gz
# Press the database
hmmpress "${interpro_hmms_dir}"/Pfam-A.hmm


# Install Java
mkdir java
curl https://launchpadlibrarian.net/172182219/openjdk-7_7u51-2.4.6.orig.tar.gz \
  --output ~/java/java.tar.gz
tar xvf ~/java/java.tar.gz -C java/
rm ~/java/java.tar.gz
cd java
# Download InterProScan


curl http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.73-104.0/interproscan-5.73-104.0-64-bit.tar.gz \
 --output ~/interproscan-5.73.tar.gz

tar xvf ~/interproscan-5.73.tar.gz
rm ~/interproscan-5.73.tar.gz
# Transdecoder
# Predicting coding regions from a transcript fasta file
# Step 1: extract the long open reading frames
TransDecoder.LongOrfs \
 -t data/ncbi_datasets/DMR_v1.0_HiC_ncbi_dataset/ncbi_dataset/data/GCF_012274545.1/rna.fna  \
 --output_dir "${transdecoder_output_dir}"/"${transdecoder_date_time}" \
  2>&1 | tee jobreports/"${transdecoder_date_time}"_DMR_rna_transdecoder_LongOrfs_report.txt

# # Step 3: predict the likely coding regions
TransDecoder.Predict \x
 -t data/ncbi_datasets/DMR_v1.0_HiC_ncbi_dataset/ncbi_dataset/data/GCF_012274545.1/rna.fna  \
 --output_dir "${transdecoder_output_dir}"/"${transdecoder_date_time}" \
  2>&1 | tee jobreports/"${transdecoder_date_time}"_DMR_rna_transdecoder_Predict_report.txt

# HMMer
# Use HMMER or InterProScan to scan proteins for domains characteristic of NF-κB components:
#     Rel Homology Domain (RHD) for NF-κB subunits.
#     ANK repeats for IκB inhibitors.
#     Kinase domains for IKK complex proteins.
# Output: Compile a list of transcripts encoding proteins with these domains.
conda activate nfkappab_evo-tools

for fasta_file in "${interpro_fasta_path}"/*.fasta
do 
 SAMPLE=$(echo "${fasta_file}" | sed "s/\.fasta//")
 domain_id=$(basename "$SAMPLE")
 echo "${domain_id}"
 clustalo -i "${fasta_file}" -o "${alignments_dir}"/"${domain_id}".afa ;
 hmmbuild --cpu "${nthreads}" \
  "${hmmprofiles_dir}"/"${domain_id}"_profile.hmm \
  "${alignments_dir}"/"${domain_id}".afa > \
  "${hmmstdout_dir}"/"${hmmer_date_time}"_"${domain_id}"_hmmbuild_stdout.txt;
 hmmsearch -A "${hmmalignments_dir}"/"${hmmer_date_time}"_"${domain_id}"_DMR_hits.sto \
  -o "${hmmstdout_dir}"/"${hmmer_date_time}"_"${domain_id}"_DMR_hmmsearch_stdout.txt \
  --domtblout "${hmmtables_dir}"/"${hmmer_date_time}"_"${domain_id}"_DMR_hmmer_dom-table.txt \
  --tblout "${hmmtables_dir}"/"${hmmer_date_time}"_"${domain_id}"_DMR_hmmer_seq-table.txt \
  --cpu "${nthreads}" \
  "${hmmprofiles_dir}"/"${domain_id}"_profile.hmm \
  "${transdecoder_output_dir}"/"${transdecoder_date_time}"/rna.fna.transdecoder.pep; 
 # Convert hits from stockholm into fasta file
 esl-reformat fasta "${hmmalignments_dir}"/"${hmmer_date_time}"_"${domain_id}"_DMR_hits.sto > \
  "${hmmfasta_dir}"/"${hmmer_date_time}"_"${domain_id}"_DMR_hits.fa;
 # Extract full sequences from domain hits
 grep -v "^#" "${hmmtables_dir}"/"${hmmer_date_time}"_"${domain_id}"_DMR_hmmer_seq-table.txt | \
  awk '{print $1}' | \
  esl-sfetch -f "${transdecoder_output_dir}"/"${transdecoder_date_time}"/rna.fna.transdecoder.pep - > \
  "${hmmfasta_dir}"/"${hmmer_date_time}"_"${domain_id}"_DMR_hits_full_seqs.fa;

done

# --domtblout <f>  : save parseable table of per-domain hits to file <f>
# --tblout <f>     : save parseable table of per-sequence hits to file <f>
# -A <f>           : save multiple alignment of all hits to file <f>

esl-sfetch --index \
 "${transdecoder_output_dir}"/"${transdecoder_date_time}"/rna.fna.transdecoder.pep

# Combine all HMM profiles
cat "${hmmprofiles_dir}"/*_profile.hmm > "${hmmprofiles_dir}"/all_custom_domain_profiles.hmm
# Press the database
hmmpress "${hmmprofiles_dir}"/all_custom_domain_profiles.hmm
# Annotate full-length reference proteins using custom HMMs
hmmscan --domtblout "${hmmtables_dir}"/"${hmmer_date_time}"_nfkb_pw_ref_all.domtblout \
 "${hmmprofiles_dir}"/all_custom_domain_profiles.hmm \
 "${full_fasta_path}"/all_selected_aa.fasta > \
 "${hmmscan_out_dir}"/nfkb_pw_ref_all.hmmscan.out

# Use a pre-built HMM to compare the results. Can it recognise the correct number of ANK repeats?
hmmscan --domtblout "${hmmtables_dir}"/"${hmmer_date_time}"_nfkb_pw_ref_all_pfam_a_prebuilt.domtblout \
 "${interpro_hmms_dir}"/Pfam-A.hmm \
 "${full_fasta_path}"/all_selected_aa.fasta > \
 "${hmmscan_out_dir}"/nfkb_pw_ref_all_pfam_a_prebuilt.hmmscan.out


hmmsearch -A "${hmmalignments_dir}"/"${hmmer_date_time}"_nfkb1_DMR_hits.sto \
  -o "${hmmstdout_dir}"/"${hmmer_date_time}"_nfkb1_DMR_hmmsearch_stdout.txt \
  --domtblout "${hmmtables_dir}"/"${hmmer_date_time}"_nfkb1_DMR_hmmer_dom-table.txt \
  --tblout "${hmmtables_dir}"/"${hmmer_date_time}"_nfkb1_DMR_hmmer_seq-table.txt \
  --cpu "${nthreads}" \
  "${interpro_hmms_dir}"/Pfam-A.hmm \
  output/fasta/dmr-nfkb1.fa

# Filter the output by e-value and domain coverage
evalue_cutoff="0.0001"
coverage_cutoff=0.5

grep -v '^#' "${hmmtables_dir}"/"${hmmer_date_time}"_nfkb_pw_ref_all.domtblout | \
 awk -v ecut="${evalue_cutoff}" -v ccut="${coverage_cutoff}" '
  {
    domain = $1
    protein = $4
    hmm_len = $3
    ali_start = $18
    ali_end = $19
    evalue = $13

    # compute domain alignment length and coverage
    align_len = (ali_end > ali_start ? ali_end - ali_start + 1 : ali_start - ali_end + 1)
    coverage = align_len / hmm_len

    if (evalue <= ecut && coverage >= ccut) {
        printf("%s\t%s\t%.2e\t%.2f\t%d\t%d\n", protein, domain, evalue, coverage, ali_start, ali_end)
    }
  }' | \
 sort -k1,1 -k5,5n  | \
 awk 'BEGIN {printf("%s\t%s\t%s\t%s\t%s\t%s", 
                "protein", "domain", "evalue", "coverage", "ali_start", "ali_end");
          print "";}1 ' > "${hmmtables_dir}"/"${hmmer_date_time}"_nfkb_pw_ref_filtered_hits.tsv






EXPECTED_DOMAINS=("Rel" "IPT" "ANK")
EXPECTED_DOMAINS=("IPR032397" "IPR002909" "IPR002110")

awk  'BEGIN {FS="\t"} NR>1{print $1"\t"$2}' \
  "${hmmtables_dir}"/"${hmmer_date_time}"_nfkb_pw_ref_filtered_hits.tsv | \
 sort -u | \
awk -v domains="$(IFS=,; echo "${EXPECTED_DOMAINS[*]}")" '
 BEGIN {
  n = split(domains, expected, ",")
 }
 {
  dom_count[$1][$2]++
 }
 END {
    for (prot in dom_count) {
      match = 1
      for (i = 1; i <= n; i++) {
        d = expected[i]
        if (dom_count[prot][d] != 1) {
          match = 0
          break
        }
      }
      status = match ? "MATCH" : "MISMATCH"
      printf("%s\t%s\t", prot, status)
      for (i = 1; i <= n; i++) {
        d = expected[i]
        count = dom_count[prot][d] ? dom_count[prot][d] : 0
        printf("%s:%d ", d, count)
      }
      print ""
    }
  }'


awk 'NR>1{print $1"\t"$2"\t"$5}' "${hmmtables_dir}"/"${hmmer_date_time}"_nfkb_pw_ref_filtered_hits.tsv | sort -k1,1 -k3,3n | \
awk '
{
  if ($1 != prev) {
    if (prev != "") print prev "\t" arch
    arch = $2
    prev = $1
  } else {
    arch = arch "," $2
  }
}
END {
  print prev "\t" arch
}'
