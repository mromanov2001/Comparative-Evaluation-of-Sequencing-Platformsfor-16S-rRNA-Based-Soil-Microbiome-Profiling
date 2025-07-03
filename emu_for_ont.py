import pandas as pd

work_dir = '/work_dir'

parse_emu_script = '/path/to/parse_emu_output.R'
parse_nanostat_script = '/path/to/parse_nanostat_output.R'

metadata = pd.read_csv(work_dir + "/metadata.csv", sep = ",", date_format = "str")
samples = metadata['barcode'].tolist()
print(samples)

rule run_16S:
    input:
        expand("{work_dir}/emu/{sample}/{sample}_rel-abundance.tsv", sample = samples, work_dir = work_dir),
        expand("{work_dir}/stat/raw/{sample}.tsv", sample = samples, work_dir = work_dir),
        expand("{work_dir}/stat/raw__trimmed__filtered/{sample}.tsv", sample = samples, work_dir = work_dir),
        expand("{work_dir}/reports/emu.tsv", work_dir = work_dir), expand("{work_dir}/reports/statistic.tsv", work_dir = work_dir)

rule get_conc_emu_table:
    input: 
        metadata = "{work_dir}/metadata.csv",
        workdir = "{work_dir}",
        emu_tables = expand("{work_dir}/emu/{sample}/{sample}_rel-abundance.tsv", sample = samples, work_dir = work_dir)
    conda: "/path/to/r_parse_libs.yaml"
    output: 
        emu = "{work_dir}/reports/emu.tsv",
        taxonomy = "{work_dir}/reports/taxonomy.tsv",
        emu_report = "{work_dir}/reports/emu_report.tsv",
    shell: ("Rscript {parse_emu_script} '{input.workdir}' '{input.metadata}' '{output.emu}' '{output.taxonomy}' '{output.emu_report}' ")

rule get_conc_stat_table:
    input:
        metadata = "{work_dir}/metadata.csv",
        workdir = "{work_dir}",
        raw = expand("{work_dir}/stat/raw/{sample}.tsv", sample = samples, work_dir = work_dir),
        filtered = expand("{work_dir}/stat/raw__trimmed__filtered/{sample}.tsv", sample = samples, work_dir = work_dir)
    conda: "/path/to/r_parse_libs.yaml"
    output: 
        "{work_dir}/reports/statistic.tsv",
    shell: ("Rscript {parse_nanostat_script} '{input.workdir}' '{input.metadata}' '{output}' ")

rule nanostat_filtered:
    input:
        raw = "{work_dir}/reads/raw__trimmed__filtered/{sample}.fastq.gz"
    threads: 10
    conda:
        "/path/to/nanostat.yaml"
    output:
        "{work_dir}/stat/raw__trimmed__filtered/{sample}.tsv"
    shell:
        "NanoStat --fastq {input} -n {output} -t {threads} --tsv"

rule nanostat_raw:
    input:
        raw = "{work_dir}/reads/raw/{sample}.fastq.gz"
    threads: 10
    conda:
        "/data15/bio/runs-maximor/yml/nanostat.yaml"
    output:
        "{work_dir}/stat/raw/{sample}.tsv"
    shell:
        "NanoStat --fastq {input} -n {output} -t {threads} --tsv"

rule emu:
    input:
        reads = "{work_dir}/reads/raw__trimmed__filtered/{sample}.fastq.gz",
        db = "/data11/bio/databases/emu_database/"
    threads: 10
    conda:
        "/path/to/emu.yaml"
    output:
        report = "{work_dir}/emu/{sample}/{sample}_rel-abundance.tsv"
    shell: "emu abundance {input.reads} --type map-ont --db {input.db} --output-dir {work_dir}/emu/{wildcards.sample} --output-basename {wildcards.sample} --keep-files --keep-counts --keep-read-assignments --threads {threads}"

rule chopper:
    input:
        reads = "{work_dir}/reads/raw__trimmed/{sample}.fastq.gz"
    threads: 10
    conda:
        "/path/to/chopper.yaml"
    output:
        temp("{work_dir}/reads/raw__trimmed__filtered/{sample}.fastq.gz")
    shell:
        "gunzip -c {input.reads} | chopper -q 10 -l 100 --maxlength 1800 | gzip > {output}"     

rule porechop:
    input:
        reads = "{work_dir}/reads/raw/{sample}.fastq.gz"
    threads: 10
    conda:
        "/path/to/porechop.yaml"
    output:
        temp("{work_dir}/reads/raw__trimmed/{sample}.fastq.gz")
    shell:
        "porechop -i {input.reads} -o {output} -t {threads}"

rule conc_files:
    input: 
        "{work_dir}/fastq_pass/{sample}"
    output: 
        temp("{work_dir}/reads/raw/{sample}.fastq.gz"),
    shell: 
        "cat {work_dir}/fastq_pass/{wildcards.sample}/*gz > {output}"
