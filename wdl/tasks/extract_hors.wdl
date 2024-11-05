version 1.0

workflow extract_hors {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "extracts HOR arrays for each chromosome in input assembly"
    }

    input {
        File assemblyFasta
        File asmToRefPaf

        File AsHorSFBedFile
        File CenSatBedFile

        String sampleName # required to be in the format HG01530_hap2

        String dockerImage="miramastoras/centromere_scripts:v0.1"
    }

    ## use alignment to CHM13 to assign HORs to chromosomes
    ## filters out discontiguous HORs, cross-chromosomal groups where not all fully are resolved
    call locate_hors {
        input:
            CenSatBedFile=CenSatBedFile,
            asmToRefPaf=asmToRefPaf,
            assemblyFasta=assemblyFasta,
            AsHorSFBedFile=AsHorSFBedFile,
            sampleID=sampleName,
            dockerImage=dockerImage
    }
    # for each chromosome, extract the hor sequence into a separate fasta
    # rename fasta header for input to centrolign
    call extract_hor_sequence {
        input:
            horArrayBed=locate_hors.horArrayBed,
            asmToRefPaf=asmToRefPaf,
            assemblyFasta=assemblyFasta,
            sampleID=sampleName,
            dockerImage=dockerImage
    }
    output {
        File horArrayBed=locate_hors.horArrayBed
        Array[File] horFastas=extract_hor_sequence.horFastas
    }
}

task locate_hors {
    input {
      File CenSatBedFile
      File asmToRefPaf
      File assemblyFasta
      File AsHorSFBedFile

      String sampleID
      String dockerImage="miramastoras/centromere_scripts:v0.1"

      Int memSizeGB=16
      Int threadCount=8
      Int diskSizeGB=64
    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        ASM_FILENAME=$(basename -- "~{assemblyFasta}")

        if [[ $ASM_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} ~{sampleID}.fasta.gz
            gunzip -f ~{sampleID}.fasta.gz
        else
            cp ~{assemblyFasta} ~{sampleID}.fasta
        fi

        samtools faidx ~{sampleID}.fasta

        python3 /opt/centromere-scripts/benchmarking/locate_hors_from_censat.py \
            -c ~{CenSatBedFile} \
            -a ~{AsHorSFBedFile} \
            -p ~{asmToRefPaf} \
            -f ~{sampleID}.fasta \
            > ~{sampleID}_hor_arrays.bed
    >>>
    output {
        File horArrayBed=glob("*hor_arrays.bed")[0]

    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task extract_hor_sequence {
    input {
      File horArrayBed
      File asmToRefPaf
      File assemblyFasta

      String sampleID
      String dockerImage="miramastoras/centromere_scripts:v0.1"

      Int memSizeGB=16
      Int threadCount=8
      Int diskSizeGB=64
    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        # get sample id without haplotype label
        SAMPLE=$(echo "~{sampleID}" | sed 's/_hap[12]//')

        # guide tree will have "hap1" and "hap2" labelled the reverse as the output from HPRC censat annotation
        if [[ "~{sampleID}" == *"hap1"* ]]; then
          PARNUM=2
          HPRC_PARENT=hap1
        else
          PARNUM=1
          HPRC_PARENT=hap2
        fi

        mkdir -p ./~{sampleID}_hor_fastas/

        for CHR in {1..22} X Y M; do
            echo "chr${CHR}"

            REGIONFILE=~{sampleID}.chr${CHR}.hor.txt
            echo $REGIONFILE
            grep -w chr${CHR} ~{horArrayBed} | awk '{ printf "%s:%d-%d\n", $1, $2+1, $3 }' > "$REGIONFILE"
            ls $REGIONFILE
            STRAND=$(grep -w chr${CHR} ~{horArrayBed} | cut -f 6)
            echo $STRAND
            if [ -s $REGIONFILE ];
            then
                echo "chr${CHR} exists, $REGIONFILE"
                # extract and add the sample name as the sequence name
                echo "extract region" `cat $REGIONFILE`
                echo "strand" $STRAND

                HORFASTA=~{sampleID}_hor_fastas/~{sampleID}.${SAMPLE}.${PARNUM}.chr${CHR}_hor_array.fasta

                samtools faidx -r $REGIONFILE ~{assemblyFasta} | sed "s/>/>$SAMPLE.$PARNUM /g" > $HORFASTA

                if [ $STRAND = "-" ];
                then
                    echo "reverse complementing sequence"
                    TMPFILE=$(mktemp)
                    /opt/centromere-scripts/data_processing_utils/fasta_to_rev_comp.py $HORFASTA > $TMPFILE
                    mv $TMPFILE $HORFASTA
                fi
            else
              echo "~{sampleID} chr${CHR} was filtered out"
            fi
        done


    >>>
    output {
        Array[File] horFastas=glob("~{sampleID}_hor_fastas/*.fasta")
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
