## Initial run of censat annotations
## this run was problematic because the repeat masker annotation is messed 
## up, but it will work for pulling active arrays (since that doesn't rely
## upon RM outputs).

cd /private/groups/hprc/qc/

mkdir -p batch1/censat
cd batch1/censat

###############################################################################
##                             Run Censat WDL                                ##
###############################################################################

wget https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/assemblies_pre_release_v0.1.index.csv


cat<<EOF>rewrite_input.py
#!/usr/bin/env python3
import csv
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, 'r', encoding='utf-8-sig') as infile, open(output_file, 'w', newline='', encoding='utf-8') as outfile:
    reader = csv.DictReader(infile)
    writer = csv.writer(outfile)
    
    # Write header
    writer.writerow(['sample_id', 'asm'])
    
    for row in reader:
        sample_id = row['sample_id']
        hap1_fa_gz = row['hap1_fa_gz']
        hap2_fa_gz = row['hap2_fa_gz']
        
        # Write hap1 row
        writer.writerow([f"{sample_id}_hap1", hap1_fa_gz])
        
        # Write hap2 row
        writer.writerow([f"{sample_id}_hap2", hap2_fa_gz])

print(f"Processed data written to {output_file}")
EOF

python3 rewrite_input.py \
    assemblies_pre_release_v0.1.index.csv \
    batch1_censat.csv

cat <<EOF > censat_input_mapping.csv 
input,type,value
centromereAnnotation.fasta,scalar,\$input.asm
centromereAnnotation.RM2Bed,scalar,/private/home/juklucas/github/alphaAnnotation/utilities/RM2Bed.py
centromereAnnotation.rDNAhmm_profile,scalar,/private/home/juklucas/github/alphaAnnotation/utilities/rDNA1.0.hmm
centromereAnnotation.AS_hmm_profile,scalar,/private/home/juklucas/github/alphaAnnotation/utilities/AS-HORs-hmmer3.4-071024.hmm
centromereAnnotation.AS_hmm_profile_SF,scalar,/private/home/juklucas/github/alphaAnnotation/utilities/AS-SFs-hmmer3.0.290621.hmm
EOF


## check that github repo is up to date
git -C /private/home/juklucas/github/alphaAnnotation pull

## check which commit we are on:
git -C /private/home/juklucas/github/alphaAnnotation log -1
# commit f554182fb3fb37d7bbb5a653a181b9338f51267e (HEAD -> main, origin/main, origin/HEAD)

mkdir input_jsons
cd input_jsons

## create input jsons for assembly cleanup
python3 /private/groups/hprc/hprc_intermediate_assembly/hpc/launch_from_table.py \
     --data_table ../batch1_censat.csv \
     --field_mapping ../censat_input_mapping.csv \
     --workflow_name censat

cd ..

mkdir -p slurm_logs

sbatch \
     --job-name=censat \
     --array=[1-312]%25 \
     --cpus-per-task=64 \
     --mem=400gb \
     --partition=long \
     /private/groups/hprc/hprc_intermediate_assembly/hpc/toil_sbatch_single_machine.sh \
     --wdl /private/home/juklucas/github/alphaAnnotation/cenSatAnnotation/centromereAnnotation.wdl \
     --sample_csv batch1_censat.csv \
     --input_json_path '../input_jsons/${SAMPLE_ID}_censat.json'


## collect results into data table    
python3 /private/groups/hprc/hprc_intermediate_assembly/hpc/update_table_with_outputs.py \
    --input_data_table batch1_censat.csv \
    --output_data_table batch1_censat_outputs.csv  \
    --json_location '{sample_id}_centromereAnnotation_outputs.json'


awk -F',' '$3 != ""' batch1_censat_outputs.csv > batch1_censat_outputs_completed.csv

cat batch1_censat_outputs.csv | wc -l            # 313
cat batch1_censat_outputs_completed.csv | wc -l  # 190

###############################################################################
##                          Upload To HPRC Bucket                            ##
###############################################################################

head -n1 batch1_censat_outputs_completed.csv
# sample_id,asm,centromeres,cenSatStrand,cenSatAnnotations,as_sf_bed,as_hor_bed,as_strand_bed,as_hor_sf_bed,RMOut

cat <<EOF > copy_linking_map.csv
column_name,destination
centromeres,reorganized/censat_test/{sample_id}
cenSatStrand,reorganized/censat_test/{sample_id}
cenSatAnnotations,reorganized/censat_test/{sample_id}
as_sf_bed,reorganized/censat_test/{sample_id}
as_hor_bed,reorganized/censat_test/{sample_id}
as_strand_bed,reorganized/censat_test/{sample_id}
as_hor_sf_bed,reorganized/censat_test/{sample_id}
RMOut,reorganized/censat_test/{sample_id}
EOF


python3 /private/groups/hprc/hprc_intermediate_assembly/hpc/misc/link_to_subfolder.py \
     --csv_file batch1_censat_outputs_completed.csv \
     --mapping_csv copy_linking_map.csv

## copy from soft linked reorg folder to new folder in Mira's directory
cp -Lr reorganized/* /private/groups/patenlab/mira/centrolign/annotations/

## update table to reflect new locations 
cat<<EOF > update_paths.py 
import sys
import re

def convert_path(old_path):
    # Skip conversion if the path doesn't match our expected pattern
    if not old_path.startswith('/private/groups/hprc/'):
        return old_path
        
    # Extract the relevant parts using regex
    pattern = r'/private/groups/hprc/qc/batch1/censat/([^/]+)/analysis/centromereAnnotation_outputs/[^/]+/(.+\.bed)'
    match = re.match(pattern, old_path)
    
    if match:
        sample_dir = match.group(1)  # e.g., HG00408_hap1
        filename = match.group(2)     # e.g., HG00408_pat_hprc_r2_v1.active.centromeres.bed
        
        # Construct new path
        return f'/private/groups/patenlab/mira/centrolign/annotations/censat_test/{sample_dir}/{filename}'
    return old_path

def main():
    # Read and output the header line unchanged
    header = sys.stdin.readline().strip()
    print(header)
    
    # Process each data line
    for line in sys.stdin:
        fields = line.strip().split(',')
        
        # Keep first two columns unchanged, convert paths in remaining columns
        new_fields = (
            fields[:2] +  # Keep sample_id and asm columns unchanged
            [convert_path(field) for field in fields[2:]]  # Convert paths in remaining columns
        )
        
        # Print the converted line
        print(','.join(new_fields))

if __name__ == '__main__':
    main()
EOF


cat batch1_censat_outputs_completed.csv \
    | python update_paths.py \
    > /private/groups/patenlab/mira/centrolign/annotations/censat_test/batch1_censat_initial_test.csv

###############################################################################
##                                   DONE                                    ##
###############################################################################