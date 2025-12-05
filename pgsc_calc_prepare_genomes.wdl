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
    }

    if (merge_chroms) {
        call merge_files {
            input:
                pgen = prepare_genomes.pgen,
                pvar = prepare_genomes.pvar,
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
        plink2 ~{prefix} ~{vcf}  \
            --allow-extra-chr \
            --chr 1-22, X, Y, XY \
            --set-all-var-ids @:#:\$r:\$a \
            ~{true="--snps-only 'just-acgt' " false="" snps_only} \
            --make-pgen multiallelics=- --out ~{basename}
    >>>

    output {
        File pgen = "~{basename}.pgen"
        File pvar = "~{basename}.pvar"
        File psam = "~{basename}.psam"
    }

    runtime {
        docker: "uwgac/pgsc_calc:0.1.0"
        #docker: "us-docker.pkg.dev/primed-cc/pgsc-calc/pgsc_calc:0.1.0"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb}G"
        cpu: "~{cpu}"
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
