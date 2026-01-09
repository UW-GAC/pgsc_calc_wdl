version 1.0

workflow pgsc_calc_prepare_genomes {
    input {
        Array[File] vcf
        Boolean merge_chroms = true
        Boolean snps_only = true
    }

    scatter (file in vcf) {
        call prepare_genomes {
            input:
                vcf = file,
                snps_only = snps_only
        }
        # Normalize alleles in the single-sample bcf
        # This changes an allele like CA/TA to C/T.
        call bcftools_norm {
            input:
                bcf = prepare_genomes.single_sample_bcf
        }
        # Convert the normalized bcf back to plink format and get the pvar file.
        call normed_pvar {
            input:
                normed_bcf = bcftools_norm.normed_bcf
        }
    }

    if (merge_chroms) {
        call merge_files {
            input:
                pgen = prepare_genomes.pgen,
                # Use the pvar file with normalized alleles.
                pvar = normed_pvar.pvar,
                psam = prepare_genomes.psam
        }
    }

    output {
        Array[File] pgen = select_first([merge_files.out_pgen, prepare_genomes.pgen])
        Array[File] pvar = select_first([merge_files.out_pvar, prepare_genomes.pvar])
        Array[File] psam = select_first([merge_files.out_psam, prepare_genomes.psam])
    }

     meta {
          author: "Stephanie Gogarten"
          email: "sdmorris@uw.edu"
    }
}


task prepare_genomes {
    input {
        File vcf
        Boolean snps_only = true
        Int mem_gb = 16
        Int cpu = 2
    }

    Int disk_size = ceil(2.5*(size(vcf, "GB"))) + 5
    String filename = basename(vcf)
    String basename = sub(filename, "[[:punct:]][bv]cf.*z?$", "")
    String prefix = if (sub(filename, ".bcf", "") != filename) then "--bcf" else "--vcf"

    command <<<
        set -e
        plink2 ~{prefix} ~{vcf}  \
            --allow-extra-chr \
            --chr 1-22, X, Y, XY \
            --set-all-var-ids @:#:\$r:\$a \
            ~{true="--snps-only 'just-acgt' " false="" snps_only} \
            --make-pgen multiallelics=- --out ~{basename}

        # Make a single sample plink for normalizing alleles in bcftools
        head -n 2 ~{basename}.psam > keep.txt
        plink2 \
            --pfile multiallelic \
            --keep keep.txt \
            --export bcf \
            --out single_sample

        plink2
    >>>

    output {
        File pgen = "~{basename}.pgen"
        File pvar = "~{basename}.pvar"
        File psam = "~{basename}.psam"
        File single_sample_bcf = "single_sample.bcf"
    }

    runtime {
        docker: "uwgac/pgsc_calc:0.1.0"
        #docker: "us-docker.pkg.dev/primed-cc/pgsc-calc/pgsc_calc:0.1.0"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb}G"
        cpu: "~{cpu}"
    }
}

task bcftools_norm {
    input {
        File bcf
    }

    String filename = basename(bcf)
    String basename = sub(filename, "[[:punct:]][bv]cf.*z?$", "")

    command <<<
        set -e
        bcftools \
            norm -a \
            ~{bcf} \
            -o ~{basename}_normed.bcf
    >>>

    output {
        File normed_bcf = "{basename}_normed.bcf"
    }

    runtime {
        docker: "SOME_DOCKER_WITH_BCFTOOLS"
    }

}

task normed_pvar {
    input {
        File normed_bcf
    }

    command <<<
        set -e
        plink2 \
            --bcf ~{normed_bcf} \
            --make-pgen \
            --out tmp
    >>>

    output {
        File new_pvar = "tmp.pvar"
    }

    runtime {
        docker: "SOME_DOCKER_WITH_BCFTOOLS_AND_PGTOOLS"
    }
}

task merge_files {
    input {
        Array[File] pgen
        Array[File] pvar
        Array[File] psam
        Int mem_gb = 16
    }

    Int disk_size = ceil(3*(size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB"))) + 10

    command <<<
        set -e -o pipefail
        cat ~{write_lines(pgen)} | sed 's/.pgen//' > pfile.txt
        plink2 --pmerge-list pfile.txt pfile \
            --merge-max-allele-ct 2 \
            --out merged
    >>>

    output {
        Array[File] out_pgen = ["merged.pgen"]
        Array[File] out_pvar = ["merged.pvar"]
        Array[File] out_psam = ["merged.psam"]
    }

    runtime {
        docker: "quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0"
        disks: "local-disk " + disk_size + " SSD"
        memory: mem_gb + " GB"
    }
}
