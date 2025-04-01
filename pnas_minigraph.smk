GENOMES, = glob_wildcards("genomes/{genome}_T2T_Chr01.fasta")

print("Total genome number ",len(GENOMES))


rule all:
    input:
        expand("graph/{asb}_graph.gfa",asb=["xiaoming"]),
        expand("remap/{asb}/{genome}_{asb}.gaf",asb=["xiaoming"],genome=GENOMES),
        expand("remap/{asb}_edge_use.tsv",asb=["xiaoming"]),
        expand("analysis/colour_node/{asb}_nodecol.tsv",asb=["xiaoming"]),
        expand("analysis/core_nonref/{asb}_core_analysis.tsv",asb=["xiaoming"])


rule construct_graph:
    input:
        ancient(expand("genomes/{genome}_T2T_Chr01.fasta",genome=GENOMES))
    output:
        "graph/{asb}_graph.gfa",
        "graph/{asb}_graph_length.tsv",
        "graph/{asb}_graph_link.tsv"
    threads: 16
    resources:
        mem_mb = 12000,
        walltime = "02:00"
    shell:
        """

        minigraph --inv no -xggs -t {threads} {input}  > {output[0]}

        awk '$1~/S/ {{ split($5,chr,":"); split($6,pos,":"); split($7,arr,":");
            print $2,length($3),chr[3],pos[3],arr[3] }}' {output[0]} > {output[1]}

        awk '$1 == "L"' {output[0]} > {output[2]}

        """

rule remap_graph:
    input:
        rules.construct_graph.output[0],
        "genomes/{genome}_T2T_Chr01.fasta"
    output:
        "remap/{asb}/{genome}_{asb}.gaf"
    threads: 2
    resources:
        mem_mb = 10000,
        walltime = "04:00"
    shell:
        """
         minigraph -t {threads} --cov -x asm {input[0]} {input[1]} > {output}
        """

def get_assemb(assemb):
    return GENOMES

rule comb_coverage:
    input:
        expand(rules.remap_graph.output,asb=["xiaoming"],genome=GENOMES)
    output:
        "remap/{asb}_coverage.tsv",
        "remap/{asb}_edge_use.tsv"
    params:
        anims = lambda wildcards: get_assemb(wildcards.asb)
    threads: 5
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """
        python scripts/comb_coverage.py -g {wildcards.asb} -a {params.anims}
        """


rule colour_node:
    input:
        rules.construct_graph.output[1],
        ancient(rules.comb_coverage.output[0])
    output:
        "analysis/colour_node/{asb}_nodecol.tsv",
        "analysis/colour_node/{asb}_nodemat.tsv"
    threads: 5,
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    params:
        assemb = lambda wildcards: get_assemb(wildcards.asb)
    shell:
        """
        Rscript scripts/colour_node.R {wildcards.asb} {params.assemb}
        """

rule identify_core_nonref:
    input:
        rules.construct_graph.output[1],
        rules.colour_node.output[1]
    output:
        "analysis/core_nonref/{asb}_core_analysis.tsv",
        "analysis/core_nonref/{asb}_nonref_analysis.tsv",
        multiext("analysis/core_nonref/{asb}_core_flex_sim", ".tsv", ".png", ".pdf"),
        multiext("analysis/core_nonref/{asb}_nonref_shared_count", ".png", ".pdf"),
        multiext("analysis/core_nonref/{asb}_nonref_shared_len", ".png", ".pdf")
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "01:00"
    script:
        "scripts/run_core_nonref.R"