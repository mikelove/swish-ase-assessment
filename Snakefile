RUNS, = glob_wildcards("/proj/milovelab/wu/bulk-ase/ase-sim/{run}_1.shuffled.fa")

SALMON = "/proj/milovelab/bin/salmon-1.4.0_linux_x86_64/bin/salmon"

rule all:
  input: expand("quants/{run}/quant.sf", run=RUNS)

rule salmon_quant:
    input:
        r1 = "/proj/milovelab/wu/bulk-ase/ase-sim/{sample}_1.shuffled.fa",
        r2 = "/proj/milovelab/wu/bulk-ase/ase-sim/{sample}_2.shuffled.fa",
        index = "bulk-ase-sim_1.4.0"
    output:
        "quants/{sample}/quant.sf"
    params:
        dir = "quants/{sample}"
    shell:
        "{SALMON} quant -i {input.index} -l A -p 12 "
        "--numGibbsSamples 20 --thinningFactor 100 "
        "-o {params.dir} -1 {input.r1} -2 {input.r2}"
