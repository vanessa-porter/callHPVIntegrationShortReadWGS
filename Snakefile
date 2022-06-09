## Load config values
configfile: "config/samples.yaml"
configfile: "config/parameters.yaml"

# path to the reference genome fasta
genome_path = config["genome_path"]

# samples
samples_dict = config["samples"]
sample_ids = samples_dict.keys()

### -------------------------------------------------------------------
### Target rule
### -------------------------------------------------------------------
rule all:
    input:
	    expand("output/{sample}/integration_sites.bed", sample=sample_ids)

### -------------------------------------------------------------------
### Rules
### -------------------------------------------------------------------

rule run_manta:
    input:
        tumour = lambda w: config["samples"][w.sample]["tumour"],
        ref = genome_path
    output:
        "output/{sample}/manta/results/variants/candidateSV.vcf.gz"
    conda: "config/conda.yaml"
    threads: 24
    shell:
        """
        configManta.py --tumorBam={input.tumour} --referenceFasta={input.ref} --runDir=output/{wildcards.sample}/manta
        output/{wildcards.sample}/manta/runWorkflow.py -j {threads}
        """

rule select_HPV:
    input:
        "output/{sample}/manta/results/variants/candidateSV.vcf.gz"
    output:
        "output/{sample}/int/hpvIntOut.vcf"
    conda: "config/conda.yaml"
    shell:
        """
        zcat {input} | grep -E "(#|HPV)" | grep -E "(#|MantaBND)" | awk 'BEGIN {{OFS=FS="\t"}} $1 !~ /^(HPV)/' > {output}
        """

## Require five reads in the PAIR_COUNT filter to count the integration site
rule read_filter:
    input:
        "output/{sample}/int/hpvIntOut.vcf"
    output:
        "output/{sample}/int/hpvIntFilt.vcf"
    conda: "config/conda.yaml"
    shell:
        """
        bcftools view -e 'PAIR_COUNT < 5' {input} > {output}
        """

rule make_bed:
    input:
        "output/{sample}/int/hpvIntFilt.vcf"
    output:
        "output/{sample}/int/hpvintSites.bed"
    conda: "config/conda.yaml"
    shell:
        """
        cat {input} | grep -v "#" | cut -f 1,2 > {output}
        """

rule site_depth:
    input:
        tumour = lambda w: config["samples"][w.sample]["tumour"],
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
        depth = "output/{sample}/int/hpvSitesDepth.txt"
    output:
        "output/{sample}/integration_sites.bed"
    conda: "config/conda.yaml"
    shell:
        """
        scripts/intSummary.R --vcf {input.vcf} --depth {input.depth} --outdir output/{wildcards.sample}
        """
