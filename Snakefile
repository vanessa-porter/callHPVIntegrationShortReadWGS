## Load config values
configfile: "config/samples.yaml"
configfile: "config/parameters.yaml"

# path to the reference genome fasta
genome_path = config["genome_path"]

# samples
samples_dict = config["bams"]
sample_ids = samples_dict.keys()

### -------------------------------------------------------------------
### Target rule
### -------------------------------------------------------------------
rule all:
    input:
	    expand("output/{sample}/.int/il_integration_events.bed", sample=sample_ids)

### -------------------------------------------------------------------
### Rules
### -------------------------------------------------------------------

rule run_manta:
    input:
        rna = lambda w: config["bams"][w.sample]["rna"],
        #normal = lambda w: config["bams"][w.sample]["normal"],
        ref = genome_path
    output:
        "output/{sample}/manta/results/variants/rnaSV.vcf.gz"
    conda: "config/conda.yaml"
    threads: 24
    shell:
        """
        configManta.py --rna --referenceFasta={input.ref} --bam {input.rna} --runDir=output/{wildcards.sample}/manta
        output/{wildcards.sample}/manta/runWorkflow.py -j {threads}
        """

rule select_HPV:
    input:
        "output/{sample}/manta/results/variants/rnaSV.vcf.gz"
    output:
        "output/{sample}/int/hpvIntOut.vcf"
    conda: "config/conda.yaml"
    shell:
        """
        zcat {input} | grep -E "(#|HPV)" | grep -E "(#|MantaBND)" | awk 'BEGIN {{OFS=FS="\t"}} $1 !~ /^(HPV)/' > {output}
        """

## Require three reads in the PAIR_COUNT filter to count the integration site
rule read_filter:
    input:
        "output/{sample}/int/hpvIntOut.vcf"
    output:
        "output/{sample}/int/hpvIntFilt.vcf"
    conda: "config/conda.yaml"
    shell:
        """
        bcftools view -e 'PAIR_COUNT < 3' {input} > {output}
        """

rule make_bed:
    input:
        "output/{sample}/int/hpvIntFilt.vcf"
    output:
        "output/{sample}/int/hpvintSites.bed"
    conda: "config/conda.yaml"
    shell:
        """
        cat {input} | grep -v "#" | cut -f 1,2 | awk '{{print $1"\t"$2"\t"$2+1}}' > {output}
        """

rule name_sites:
    input:
        "output/{sample}/int/hpvintSites.bed"
    output:
        "output/{sample}/il_integration_sites.bed"
    conda: "config/conda.yaml"
    shell:
        """
        scripts/hpvSiteName.R --bed={input} --out={output}
        """

rule make_events:
    input:
        "output/{sample}/il_integration_sites.bed"
    output:
        "output/{sample}/int/integration_events.bed"
    conda: "config/conda.yaml"
    shell:
        """
        bedtools merge -d 500000 -o collapse -c 4 -i {input} > {output}
        """

rule site_depth:
    input:
        tumour = lambda w: config["bams"][w.sample]["rna"],
        sites = "output/{sample}/int/hpvintSites.bed"
    output:
        "output/{sample}/int/hpvSitesDepth.txt"
    conda: "config/conda.yaml"
    shell:
        """
        samtools depth -a -b {input.sites} {input.tumour} > {output}
        """

rule summary:
    input:
        vcf = "output/{sample}/int/hpvIntFilt.vcf",
        depth = "output/{sample}/int/hpvSitesDepth.txt", 
        bed = "output/{sample}/int/integration_events.bed"
    output:
        "output/{sample}/il_integration_events.bed"
    conda: "config/conda.yaml"
    shell:
        """
        scripts/intSummary.R --vcf {input.vcf} --bed={input.bed} --depth {input.depth} --outdir output/{wildcards.sample}
        """

rule mv_other_files:
    input:
        "output/{sample}/il_integration_events.bed"
    output:
        "output/{sample}/.int/il_integration_events.bed"
    conda: "config/conda.yaml"
    shell:
        """
        mv output/{wildcards.sample}/int output/{wildcards.sample}/.int
        touch {output}
        """
