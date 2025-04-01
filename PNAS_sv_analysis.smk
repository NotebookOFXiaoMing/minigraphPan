ASB , = glob_wildcards("graph/{asb}_graph.gfa")
GENOMES, = glob_wildcards("genomes/{genome}_T2T_Chr01.fasta")

svlist = ["biallelic", "multiallelic"]

rule all:
    input:
        expand("analysis/bubble/{asb}_bubble.tsv",asb=ASB),
        expand("analysis/bubble/{asb}_biallelic_sv.tsv",asb=ASB),
        expand("analysis/bubble/{asb}_bialsv_seq.fa",asb=ASB),
        expand("analysis/bubble/{asb}_multiallelic_sv.tsv",asb=ASB),
        expand("analysis/bubble/{asb}_multisv_seq.fa",asb=ASB),
        expand("analysis/bubble/{asb}_path_trace.tsv",asb=ASB),
        expand("analysis/bubble/{asb}_{svtype}_sv_viz.pdf",asb=ASB,svtype=svlist),
        expand("analysis/bubble/{asb}_nonrefsv.fa",asb=ASB),
        expand("analysis/bubble/{asb}_nonrefsv_woflank.fa",asb=ASB),
        expand("extended_ref/{asb}_Jack_T2T_Chr01_plus_nonrefsv.fa",asb=ASB)


rule identify_bubble:
    input:
        "graph/{asb}_graph.gfa"
    output:
        "analysis/bubble/{asb}_bubble.tsv",
        "analysis/bubble/{asb}_biallelic_bubble.tsv",
        "analysis/bubble/{asb}_multiallelic_bubble.tsv",
        "analysis/bubble/{asb}_bubble.bed"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """

        gfatools bubble {input} > {output[0]}

        awk '$5==2 {{ print $1,$2,$4,$5,$12 }}' {output[0]} > {output[1]}

        awk '$5>2 && $5 < 8 {{ print $1,$2,$4,$5,$12 }}' {output[0]} > {output[2]}

        awk '{{ print $1,$2,$2+1,$1"_"$2 }}' OFS="\t" {output[0]} > {output[3]}
        """

rule collect_biallelic_sv:
    input:
        "graph/{asb}_graph_length.tsv",
        rules.identify_bubble.output[1]
        #"analysis/bubble/{asb}_biallelic_bubble.tsv"
    output:
        "analysis/bubble/{asb}_biallelic_sv.tsv"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """

        python scripts/get_bialsv.py -a {wildcards.asb} > {output}

        """

rule extract_bialsv:
    input:
        "graph/{asb}_graph.gfa",
        rules.collect_biallelic_sv.output
    output:
        "analysis/bubble/{asb}_bialsv_seq.fa",
        "analysis/bubble/{asb}_bialsv_stat.tsv"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "00:30"
    shell:
        """

        python scripts/get_bialseq.py -a {wildcards.asb}

        """

rule collect_multiallelic_sv:
    input:
        "graph/{asb}_graph_length.tsv",
        "analysis/bubble/{asb}_multiallelic_bubble.tsv",
        "analysis/colour_node/{asb}_nodecol.tsv"
    output:
        "analysis/bubble/{asb}_multiallelic_sv.tsv"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """
        python scripts/get_multisv.py -a {wildcards.asb} > {output}
        """

rule extract_multisv:
    input:
        "graph/{asb}_graph.gfa",
        rules.collect_multiallelic_sv.output
    output:
        "analysis/bubble/{asb}_multisv_seq.fa",
        "analysis/bubble/{asb}_multisv_stat.tsv"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "00:30"
    shell:
        """
        python scripts/get_multiseq.py -a {wildcards.asb}
        """
def get_assemb(assemb):
    return GENOMES


rule trace_paths:
    input:
        "remap/{asb}_edge_use.tsv",
        "analysis/bubble/{asb}_biallelic_sv.tsv",
        "analysis/bubble/{asb}_multiallelic_sv.tsv"
    output:
        "analysis/bubble/{asb}_path_trace.tsv"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "01:00"
    params:
        assemb = lambda wildcards: get_assemb(wildcards.asb)
    shell:
        """

        python scripts/trace_path.py -g {wildcards.asb} -a {params.assemb} > {output}

        """

rule visualize_sv:
    input:
        rules.collect_biallelic_sv.output,
        rules.collect_multiallelic_sv.output,
        #rules.construct_graph.output
    output: "analysis/bubble/{asb}_{svtype}_sv_viz.pdf"
    threads: 10
    params:
        assemb = lambda wildcards: get_assemb(wildcards.asb)
    # envmodules:
    #     "gcc/4.8.5",
    #     "graphviz/2.40.1",
    #     "python_cpu/3.7.4"
    resources:
        mem_mb = 1000,
        walltime = "01:00"
    shell:
        """

        python visualize/sv_viz.py -g {wildcards.asb} -c {params.assemb} -m {wildcards.svtype}

        """


rule combine_sv:
    input:
        #"analysis/bubble/{asb}_bialsv_seq.fa",
        #"analysis/bubble/{asb}_multisv_seq.fa"
        rules.extract_bialsv.output[0],
        rules.extract_multisv.output[0]
    output:
        "analysis/bubble/{asb}_nonrefsv.fa"
    shell:
        """
            cat {input} > {output}
        """

rule combine_sv_woflank:
    input:
        bialfile = rules.extract_bialsv.output[1],
        multifile = rules.extract_multisv.output[1],
        graphfile = "graph/{asb}_graph.gfa"
    output:
        "analysis/bubble/{asb}_nonrefsv_woflank.fa"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    shell:
        """

        python scripts/get_sv_woflank.py \
                -b {input.bialfile} \
                -m {input.multifile} \
                -g {input.graphfile} \
                -o {output}

        """


#localrules: create_full_ref
rule create_full_ref:
    input:
        genome = "genomes/Jack_T2T_Chr01.fasta",
        sv = rules.combine_sv.output
    output: "extended_ref/{asb}_Jack_T2T_Chr01_plus_nonrefsv.fa"
    shell:
        """

        cat {input.genome} {input.sv} > {output}

        """