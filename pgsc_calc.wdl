version 1.0

import "pgsc_calc_prepare_genomes.wdl" as prep

workflow pgsc_calc {
    input {
        Array[File]? vcf
        Array[File]? pgen
        Array[File]? pvar
        Array[File]? psam
        Array[String] chromosome = [""]
        String target_build = "GRCh38"
        File? pgs_id
        File? scorefile
        File? ancestry_ref_panel
        String sampleset_name = "cohort"
        Array[String]? arguments
    }

    if (defined(vcf)) {
        scatter (file in select_first([vcf, ""])) {
            call prep.prepare_genomes {
                input:
                    vcf = file
            }
        }
    }

    if (defined(scorefile)) {
        call add_scorefile_header {
            input:
                scorefile = select_first([scorefile, ""])
        }
    }

    call pgsc_calc_nextflow {
        input:
            pgen = select_first([prepare_genomes.pgen, pgen]),
            pvar = select_first([prepare_genomes.pvar, pvar]),
            psam = select_first([prepare_genomes.psam, psam]),
            chromosome = chromosome,
            scorefile = select_first([add_scorefile_header.scorefile_hdr, scorefile]),
            pgs_id = pgs_id,
            target_build = target_build,
            ancestry_ref_panel = ancestry_ref_panel,
            sampleset = sampleset_name,
            arguments = arguments
    }

    output {
        Array[File] match_files = pgsc_calc_nextflow.match_files
        Array[File] score_files = pgsc_calc_nextflow.score_files
        Array[File] log_files = pgsc_calc_nextflow.log_files
    }

     meta {
          author: "Stephanie Gogarten"
          email: "sdmorris@uw.edu"
    }
}


task add_scorefile_header {
    input {
        File scorefile
        String pgs_name = "unknown"
        String pgs_id = "unknown"
        String trait_reported = "unknown"
        String genome_build = "GRCh38"
    }

    command <<<
        R << RSCRIPT
            outfile <- "scorefile_with_header.txt"
            chk <- readLines("~{scorefile}", n=100)
            if (!any(stringr::str_detect(chk, "^#genome_build"))) {
                header <- c(
                    "#pgs_name=~{pgs_name}",
                    "#pgs_id=~{pgs_id}",
                    "#trait_reported=~{trait_reported}",
                    "#genome_build=~{genome_build}"
                )
                writeLines(header, outfile)
            }
            dat <- readr::read_tsv("~{scorefile}", comment = "#")
            readr::write_tsv(dat, outfile, append=TRUE, col_names=TRUE)
        RSCRIPT
    >>>

    output {
        File scorefile_hdr = "scorefile_with_header.txt"
    }

    runtime {
        docker: "rocker/tidyverse:4"
        disks: "local-disk 16 SSD"
        memory: "8G"
    }
}


task pgsc_calc_nextflow {
    input {
        Array[File] pgen
        Array[File] pvar
        Array[File] psam
        Array[String] chromosome
        String target_build
        String? pgs_id
        File? scorefile
        File? ancestry_ref_panel
        String sampleset
        Array[String]? arguments
        Int disk_gb = 128
        Int mem_gb = 64
        Int cpu = 16
    }


    command <<<
        set -e -o pipefail

        Rscript -e "\
        files <- readLines('~{write_lines(pgen)}'); \
        chrs <- readLines('~{write_lines(chromosome)}'); \
        stopifnot(length(files) == length(chrs)); \
        file_prefix <- sub('.pgen$', '', files); \
        sampleset <- tibble::tibble(sampleset = '~{sampleset}', path_prefix=file_prefix, chrom=chrs, format='pfile'); \
        readr::write_csv(sampleset, 'samplesheet.csv'); \
        "

        nextflow run pgscatalog/pgsc_calc -r v2.0.0-alpha.5 -profile conda \
            --input samplesheet.csv \
            --target_build ~{target_build} \
            ~{"--pgs_id " + pgs_id} \
            ~{"--scorefile " + scorefile} \
            ~{"--run_ancestry " + ancestry_ref_panel} \
            ~{sep=" " arguments}
    >>>

    output {
        File samplesheet = "samplesheet.csv"
        Array[File] match_files = glob("results/~{sampleset}/match/*")
        Array[File] score_files = glob("results/~{sampleset}/score/*")
        Array[File] log_files = glob("results/pipeline_info/*")
        File nextflow_log = ".nextflow.log"
    }

    runtime {
        docker: "uwgac/pgsc_calc:0.1.0"
        #docker: "us-docker.pkg.dev/primed-cc/pgsc-calc/pgsc_calc:0.1.0"
        disks: "local-disk ~{disk_gb} SSD"
        memory: "~{mem_gb}G"
        cpu: "~{cpu}"
    }
}
