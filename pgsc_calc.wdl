version 1.0

workflow pgsc_calc {
    input {
        
    }

    call pgsc_calc_nextflow {

    }

    output {
        
    }

     meta {
          author: "Stephanie Gogarten"
          email: "sdmorris@uw.edu"
     }
}

task pgsc_calc_nextflow {
    input {
        Int mem_gb = 16
        Int cpu = 2
    }

    command <<<
        nextflow run pgscatalog/pgsc_calc -r v2.0.0-alpha.5 -profile test,conda
    >>>

    output {

    }

    runtime {
        docker: "uwgac/pgsc_calc:0.1.0"
        memory: "~{mem_gb}G"
        cpu: "~{cpu}"
    }
}
