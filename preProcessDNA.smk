rule all:
    input:
        "cleaned_forward.fastq",
        "cleaned_reverse.fastq"

rule preprocess_reads:
    input:
        input_forward="testdata/test_forward.fastq",
        input_reverse="testdata/test_reverse.fastq",
        adapters="testdata/adapters.fna"
    output:
        cleaned_forward="cleaned_forward.fastq",
        cleaned_reverse="cleaned_reverse.fastq"
    container: "docker://ghcr.io/bwbioinfo/bbtools-docker-cwl:latest"
    shell:
        """
        bash /opt/scripts/preProcessDNA.sh {input.input_forward} {input.input_reverse} {input.adapters} {output.cleaned_forward} {output.cleaned_reverse}
        """
