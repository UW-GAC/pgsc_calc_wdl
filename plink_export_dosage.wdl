version 1.0

workflow plink_export_dosage {
    input {
        File pgen
        File psam
        File pvar
        File samples
        File variants
    }

    call export_dosage {
        input:
            pgen = pgen,
            psam = psam,
            pvar = pvar,
            samples = samples,
            variants = variants
    }

    output {
        File dosage = export_dosage.dosage
    }
}


task export_dosage {
    input {
        File pgen
        File psam
        File pvar
        File samples
        File variants
        String prefix = "out"
        Int mem_gb = 16
        Int cpu = 2
    }

    command <<<
        plink2 --pgen ~{pgen} --pvar ~{pvar} --psam ~{psam} \
            --keep samples --extract variants \
            --export Av \
            --out ~{prefix}
    >>>

    Int disk_size = ceil(2*(size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB"))) + 10

    output {
        File dosage = "out.traw"
    }

    runtime {
        docker: "quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb}G"
        cpu: "~{cpu}"
    }
}
