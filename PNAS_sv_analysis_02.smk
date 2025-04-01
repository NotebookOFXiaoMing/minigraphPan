ASB , = glob_wildcards("graph/{asb}_graph.gfa")
GENOMES, = glob_wildcards("genomes/{genome}_T2T_Chr01.fasta")

svlist = ["biallelic", "multiallelic"]


rule all:
    input:
        expand("analysis/bubble/{asb}_left_breakpoints.bed",asb=ASB),
        expand("analysis/bubble/{asb}_right_breakpoints.bed",asb=ASB),


rule create_breakpoint_bed:
    input:
        bubble_file = "analysis/bubble/{asb}_bubble.tsv",
        bialsv_file = "analysis/bubble/{asb}_biallelic_sv.tsv",
        multisv_file = "analysis/bubble/{asb}_multiallelic_sv.tsv"
    output:
        left_bed = "analysis/bubble/{asb}_left_breakpoints.bed",
        right_bed = "analysis/bubble/{asb}_right_breakpoints.bed"
    run:
        import re

        # get right breakpoint
        # my chromosome name sampleid#1#Chr01
        right_bp = {}
        with open(input["bubble_file"]) as infile:
            for line in infile:
                chromo, left_side, right_side, *_ = line.strip().split()
                # only get the numeric part of the chromo
                #chromo = [int(x) for x in chromo.split("_") if re.search(r"\d+", x)][0]
                chromo = int(re.sub("^0","",chromo.split("#")[-1][-2:]))
                #print(chromo)
                right_bp[f"{chromo}_{left_side}"] = right_side

        def wrote_sv_bed(line, left_file, right_file, mutype="biallelic"):
            if mutype == "biallelic":
                # 1 165873 AltDel 497 2 s1 s2 s133016 s3
                line_comp = line.strip().split()
                chromo, leftcoord = line_comp[:2]
                print(chromo)
                chromo = re.sub("^0","",chromo.split("#")[-1][-2:])
                print(chromo)
                start_node, stop_node = line_comp[5], line_comp[-1]
                sv_comp = f"{chromo}_{leftcoord}"
                leftcoord = int(leftcoord)
                # sv_comp, *sv_rest = line.strip().split()
                # _, *chromo, leftcoord = sv_comp.split("_")
                # only get the numeric part as the chromosome id
                # chromo = [int(x) for x in chromo if re.search(r"\d+", x)][0]
                # start_node, *_, stop_node = sv_rest[-2].split(",")
            if mutype == "multiallelic":
                # 1_535561        1898    5138    AltIns  s33,s60173,s36
                line_comp = line.strip().split()
                print(line_comp)
                chromo, leftcoord = line_comp[0].split("#")[-1].split("_")
                print(chromo)
                chromo = re.sub("Chr0{0,1}","",chromo)
                leftcoord = int(leftcoord)
                sv_comp = f"{chromo}_{leftcoord}"
                start_node, *_, stop_node = line_comp[-1].split(",")
            # write the bed file
            left_file.write(
                (f"{chromo}\t{leftcoord-1}\t{leftcoord+1}\t{start_node}\t{stop_node}\t{sv_comp}\n"))
            svid = f"{chromo}_{leftcoord}"
            rightcoord = int(right_bp[svid])
            right_file.write(
                (f"{chromo}\t{rightcoord-1}\t{rightcoord+1}\t{start_node}\t{stop_node}\t{sv_comp}\n"))

        with open(input["bialsv_file"]) as bialfile, open(input["multisv_file"]) as multifile:
            with open(output["left_bed"], "a") as left_file, open(output["right_bed"], "a") as right_file:
                for line in bialfile:
                    # process the biallelic breakpoints
                    wrote_sv_bed(line, left_file, right_file, mutype="biallelic")
            # process the multiallelic breakpoints
                sv_processed = []
                mutlist = []
                for line in multifile:
                    svid = line.strip().split()[0]
                    # sv_comp, *sv_rest = line.strip().split()
                    # _, *chromo, leftcoord = sv_comp.split("_")
                    # # only get the numeric part as the chromosome id
                    # chromo = [int(x) for x in chromo if re.search(r"\d+", x)][0]
                    # svid = f"{chromo}_{leftcoord}"
                    if svid not in sv_processed:
                        sv_processed.append(svid)
                        wrote_sv_bed(line, left_file, right_file, mutype="multiallelic")
