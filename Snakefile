# Most of the steps of this pipeline are designed to be run in the cluster
# But with some modifications could be used in other server/infraestructure

from collections import defaultdict
import os

#FOLDERS
METADATA = "data/metadata/"
RAW = "data/raw/"
PROCESS = "data/process/"
REFERENCES = "data/references/"
RESULTS = "results/"

#FOLDERS, SPECIFIC TO MY CONFIGURATION
DBS = "/hpcudd/home/jugalde/storage/databases/"
path_staph_genomes = "/hpcudd/home/jugalde/storage/mrsa_chile/genomes_from_staphophia/juan-sts/contigs/"


# Processing functions
def get_ncbi_file_by_type(input_folder, sample_list, file_extension):
    sample_files = defaultdict()

    for sample in sample_list:
        for file in os.listdir(input_folder + sample):
            if file_extension == "fna":
                if file.endswith("_genomic.fna.gz") and not file.endswith("_from_genomic.fna.gz"):
                    sample_files[sample] = file

            elif file_extension == "faa":
                if file.endswith("protein.faa.gz"):
                    sample_files[sample] = file

            elif file_extension == "gff":
                if file.endswith("_genomic.gff.gz"):
                    sample_files[sample] = file

    return sample_files

def fix_names(input_list):

    return_list = list()

    for name in input_list:
        prefix = "_".join(name.split("_")[:-1])
        return_list.append(prefix)

    return return_list

# Process the large db dump of Staph genomes

staphophia_ST_genomes = defaultdict(list)

for line in open("data/juan-sts.txt", 'r'):
    line = line.rstrip()
    if line.startswith("sample"):
        continue

    elements = line.split("\t")

    staphophia_ST_genomes[str(elements[4])].append(elements[1])


# Load the RGI database
if not os.path.exists("localDB"):
    os.system("rgi load --card_json /hpcudd/home/jugalde/storage/databases/rgi/card.json --local")

# Sample list
SAMPLES = set() # MRSA from this project
chile_samples = set() # MRSA genomes from NCBI, from Chile
reference_samples = set() # References based on paper XXX
staphopia_samples = set()

for file in os.listdir(RAW + "fastq/"):
    sample_name = file.split("_")[0]
    SAMPLES.add(sample_name)

for file in os.listdir(RAW + "external_genomes/chile_data"):
    chile_samples.add(file)

for file in os.listdir(RAW + "external_genomes/reference_genomes"):
    reference_samples.add(file)

for file in os.listdir(RAW + "staphopia_genomes"):
    sample_name = file.split(".")[0]
    staphopia_samples.add(sample_name)


chile_genomes = get_ncbi_file_by_type(RAW + "external_genomes/chile_data/", chile_samples, "fna")
reference_genomes = get_ncbi_file_by_type(RAW + "external_genomes/reference_genomes/", reference_samples, "fna")

chile_proteins = get_ncbi_file_by_type(RAW + "external_genomes/chile_data/", chile_samples, "faa")
reference_proteins = get_ncbi_file_by_type(RAW + "external_genomes/reference_genomes/", reference_samples, "faa")

chile_gff = get_ncbi_file_by_type(RAW + "external_genomes/chile_data/", chile_samples, "gff")
reference_gff = get_ncbi_file_by_type(RAW + "external_genomes/reference_genomes/", reference_samples, "gff")

# Create paired lists for the gff3 creation
chile_combinations = list()
reference_combinations = list()

for entry in chile_gff:
    gff = RAW + "external_genomes/chile_data/" + entry + "/" + chile_gff[entry]
    fna = RAW + "external_genomes/chile_data/" + entry + "/" + chile_genomes[entry]

    chile_combinations.append([gff, fna, entry])

for entry in reference_genomes:
    gff = RAW + "external_genomes/reference_genomes/" + entry + "/" + reference_gff[entry]
    fna = RAW + "external_genomes/reference_genomes/" + entry + "/" + reference_genomes[entry]

    reference_combinations.append([gff, fna, entry])

# Rename GFF files for the chile and reference data. This should be a rule, but I couldn't figure out an elegant
# way of doing this

def rename_fna_files(input_sample, fna_dictionary, folder, output_folder):
    output_file = output_folder + input_sample + ".fna.gz"
    input_file = folder + "/" + input_sample + "/" + fna_dictionary[input_sample]
    os.system("cp " + input_file + " " + output_file)
    os.system("gunzip " + output_file)

all_reference_samples = list(chile_gff.keys()) + list(reference_gff.keys())


def get_gff_study(folder, mlst_file, request_st):
    gff_list = list()

    for line in open(mlst_file, 'r'):
        line = line.rstrip()

        st_type = str(line.split("\t")[2])

        if st_type == request_st:
            sample_name = line.split("\t")[0].split("/")[-2]

            full_path = folder + sample_name + "/" + sample_name + ".gff"

            gff_list.append(full_path)

    return gff_list


# Run Shovill for each of the genomes

rule all:
    input:
        # Run BacMet
        expand(RESULTS + "BacMet/{sample}.tsv", sample=SAMPLES),
        expand(RESULTS + "BacMet_staphopia_ref/{sample}.tsv", sample=staphopia_samples),
        expand(RESULTS + "BacMet_references/{sample}.tsv", sample=reference_samples|chile_samples),

        # Run MLST
        RESULTS + "MLST_results.tsv",
        RESULTS + "Chile_MLST_results.tsv",
        RESULTS + "References_MLST_results.tsv",

        # Run RGI
        expand(RESULTS + "rgi/{sample}.txt", sample=SAMPLES),

        # Run Prokka
        expand(PROCESS + "prokka/{sample}/{sample}.gff", sample=SAMPLES),
        expand(PROCESS + "prokka_staphopia/{sample}/{sample}.gff", sample=staphopia_samples),

        # Run Roary
        #RESULTS + "roary/summary_statistics.txt",
        RESULTS + "roary_mafft/summary_statistics.txt",

        # Rorary with Reference genomes

        #RESULTS + "roary_st8/summary_statistics.txt",
        #RESULTS + "roary_st105/summary_statistics.txt",
        #RESULTS + "roary_st30/summary_statistics.txt",
        #RESULTS + "roary_st5/summary_statistics.txt",
        #RESULTS + "roary_st72/summary_statistics.txt"



    # Commented until I find a way to generate fastg files from the gfa
     #   expand(RESULTS + "recycler/{sample}.cycs.fasta", sample=SAMPLES)

rule run_shovill:
    input:
        R1 = RAW + "fastq/{sample}_1_trimmed.fastq",
        R2 = RAW + "fastq/{sample}_2_trimmed.fastq",
        R1_UN = RAW + "fastq/{sample}_U1_trimmed.fastq",
        R2_UN = RAW + "fastq/{sample}_U2_trimmed.fastq"

    output:
        contigs = PROCESS + "shovill/{sample}/contigs.fa"
    params:
        output_dir = PROCESS + "shovill/{sample}"

    shell:
        "shovill --outdir {params.output_dir} --R1 {input.R1} --R2 {input.R2} --depth 100 --gsize 3000000 --cpus 20 --force"

rule run_recycler:
    input:
        R1 = RAW + "fastq/{sample}_1_trimmed.fastq",
        R2 = RAW + "fastq/{sample}_2_trimmed.fastq",
        contigs = PROCESS + "shovill/{sample}/contigs.fa",
        graph = PROCESS + "shovill/{sample}/contigs.gfa"

    output:
        recycler_output = RESULTS + "recycler/{sample}.cycs.fasta",
        mapping_results = PROCESS + "bwa_mapping/{sample}.primary.sort.bam"

    params:
        prefix = "{sample}",
        bwa_mem_file = PROCESS + "bwa_mapping/{sample}.bam",
        sam_view = PROCESS + "bwa_mapping/{sample}.primary.bam",
        sam_sort = PROCESS + "bwa_mapping/{sample}.primary.sort.bam",
        recycler_folder = RESULTS + "recycler"

    shell:
        """
        bwa index {input.contigs}
        bwa mem {input.contigs} {input.R1} {input.R2} | samtools view -buS - > {params.bwa_mem_file}
        samttols view -bF 0x800 {params.bwa_mem_file} -> {params.sam_view}
        samtools sort {params.sam_view} -o {params.sam_sort}
        samtools index {params.sam_sort}
        recycle.py -g {input.graph} -k 127 -b {params.sam_sort} -i True -o {params.recycler_folder}
        rm -f {params.bwa_mem_file}
        rm -f {params.sam_view}
        """

rule run_prokka:
    input:
        contigs = PROCESS + "shovill/{sample}/contigs.fa"

    output:
        gff = PROCESS + "prokka/{sample}/{sample}.gff"

    params:
        output_dir = PROCESS + "prokka/{sample}",
        prefix = "{sample}"

    shell:
        "prokka --outdir {params.output_dir} --prefix {params.prefix} --locustag {params.prefix} --addgenes --mincontiglen 200 --genus Staphylococcus --species aureus --strain {params.prefix} --kingdom Bacteria --rfam --cpus 20 --force {input.contigs}"


rule mlst:
    input:
        fna_files = expand(PROCESS + "prokka/{sample}/{sample}.fna", sample=SAMPLES)

    output:
        mlst_output = RESULTS + "MLST_results.tsv"

    shell:
        """
        mlst --scheme saureus {input.fna_files} > {output.mlst_output}
        """

rule run_mlst_chile:
    input:
        fna_file = expand(RAW + "external_genomes/chile_data/{sample}/{file}", zip, sample=chile_genomes.keys(), file=chile_genomes.values())

    output:
        mlst_output = RESULTS + "Chile_MLST_results.tsv"

    shell:
        """
        mlst --scheme saureus {input.fna_file} > {output.mlst_output}
        """

rule run_mlst_references:
    input:
        fna_file = expand(RAW + "external_genomes/reference_genomes/{sample}/{file}", zip, sample=reference_genomes.keys(), file=reference_genomes.values())

    output:
        mlst_output = RESULTS + "References_MLST_results.tsv"

    shell:
        """
        mlst --scheme saureus {input.fna_file} > {output.mlst_output}
        """

rule run_bacmet:
    input:
         sample_contig = PROCESS + "prokka/{sample}/{sample}.faa"

    output:
          bacmet_results = RESULTS + "BacMet/{sample}.tsv"

    params:
        dmnd_file = DBS + "bacmet/BacMet2.dmnd"

    shell:
         """
         diamond blastp -p 24 -q {input.sample_contig} -e 0.00001 -d {params.dmnd_file} -o {output.bacmet_results} --id 90
         """

rule run_bacmet_staphopia_ref:
    input:
         sample_contig = PROCESS + "prokka_staphopia/{sample}/{sample}.faa"

    output:
          bacmet_results = RESULTS + "BacMet_staphopia_ref/{sample}.tsv"

    params:
        dmnd_file = DBS + "bacmet/BacMet2.dmnd"

    threads: 27

    shell:
         """
         diamond blastp -p {threads} -q {input.sample_contig} -e 0.00001 -d {params.dmnd_file} -o {output.bacmet_results} --id 90
         """

rule run_bacmet_references:
    input:
         sample_contig = PROCESS + "prokka_references/{sample}/{sample}.faa"

    output:
          bacmet_results = RESULTS + "BacMet_references/{sample}.tsv"

    params:
        dmnd_file = DBS + "bacmet/BacMet2.dmnd"

    threads: 27

    shell:
         """
         diamond blastp -p {threads} -q {input.sample_contig} -e 0.00001 -d {params.dmnd_file} -o {output.bacmet_results} --id 90
         """

rule run_prokka_staphopia:
    input:
        contigs = RAW + "staphopia_genomes/{sample}.fna"

    output:
        gff = PROCESS + "prokka_staphopia/{sample}/{sample}.gff"

    params:
        output_dir = PROCESS + "prokka_staphopia/{sample}",
        prefix = "{sample}"

    shell:
        "prokka --outdir {params.output_dir} --prefix {params.prefix} --locustag {params.prefix} --addgenes --mincontiglen 200 --genus Staphylococcus --species aureus --strain {params.prefix} --kingdom Bacteria --cpus 20 --force {input.contigs}"


rule copy_fna_references:
    output:
        output_files = expand(RAW + "external_genomes/fna_files/{sample}.fna", sample=all_reference_samples)
    run:
        for sample in chile_gff.keys():
            rename_fna_files(sample, chile_genomes, RAW + "external_genomes/chile_data/", RAW + "external_genomes/fna_files/")

        for sample in reference_gff.keys():
            rename_fna_files(sample, reference_genomes, RAW + "external_genomes/reference_genomes/", RAW + "external_genomes/fna_files/")

rule run_prokka_references:
    input:
        contigs = RAW + "external_genomes/fna_files/{sample}.fna"

    output:
        gff = PROCESS + "prokka_references/{sample}/{sample}.gff"

    params:
        output_dir = PROCESS + "prokka_references/{sample}",
        prefix = "{sample}"

    shell:
        "prokka --outdir {params.output_dir} --prefix {params.prefix} --locustag {params.prefix} --addgenes --mincontiglen 200 --genus Staphylococcus --species aureus --strain {params.prefix} --kingdom Bacteria --cpus 20 --force {input.contigs}"


rule run_roary:
     input:
         gff_files = expand(PROCESS + "prokka/{sample}/{sample}.gff", sample=SAMPLES) + expand(PROCESS + "prokka_staphopia/{sample}/{sample}.gff", sample=staphopia_samples) + expand(PROCESS + "prokka_references/{sample}/{sample}.gff", sample=chile_gff.keys()) + expand(PROCESS + "prokka_references/{sample}/{sample}.gff", sample=reference_gff.keys())

     output:
         roary_output = RESULTS + "roary/summary_statistics.txt"

     params:
         output_folder = RESULTS + "roary"

     threads: 27

     shell:
         """
         rm -rf {params.output_folder}
         roary -e -p {threads} -f {params.output_folder} {input.gff_files}
         """

rule run_rgi:
    input:
        faa = PROCESS + "prokka/{sample}/{sample}.faa"

    output:
        rgi_output = RESULTS + "rgi/{sample}.txt"

    wildcard_constraints:
        sample = "Sa\d+"

    params:
        output_file = RESULTS + "rgi/{sample}"

    shell:
        """
        rgi main --input_sequence {input.faa} --output_file {params.output_file} --input_type protein --local --include_loose --clean 
        """

rule annotate_staphophia_genomes:
    input:
        contigs = path_staph_genomes + "{sample}.contigs.fasta"

    output:
        gff = PROCESS + "prokka_STS_staphophia/{sample}/{sample}.gff"

    params:
        output_dir = PROCESS + "prokka_STS_staphophia/{sample}",
        prefix = "{sample}"

    conda:
        "prokka.yml"

    threads: 27

    shell:
        "prokka --outdir {params.output_dir} --prefix {params.prefix} --locustag {params.prefix} --addgenes --mincontiglen 200 --genus Staphylococcus --species aureus --strain {params.prefix} --kingdom Bacteria --cpus {threads} --force {input.contigs}"

rule run_roary_mafft:
     input:
         gff_files = expand(PROCESS + "prokka/{sample}/{sample}.gff", sample=SAMPLES) + expand(PROCESS + "prokka_staphopia/{sample}/{sample}.gff", sample=staphopia_samples) + expand(PROCESS + "prokka_references/{sample}/{sample}.gff", sample=chile_gff.keys()) + expand(PROCESS + "prokka_references/{sample}/{sample}.gff", sample=reference_gff.keys())

     output:
         roary_output = RESULTS + "roary_mafft/summary_statistics.txt"

     params:
         output_folder = RESULTS + "roary_mafft"

     threads: 27

     shell:
         """
         rm -rf {params.output_folder}
         roary -e -n -p {threads} -f {params.output_folder} {input.gff_files}
         """

rule run_roary_ST_72:
     input:
          gff_files = [PROCESS + "prokka_STS_staphophia/" + accession + "/" + accession + ".gff" for accession in staphophia_ST_genomes["72"]] + get_gff_study(PROCESS + "prokka/", RESULTS + "MLST_results.tsv", "72")

     output:
         roary_output = RESULTS + "roary_st72/summary_statistics.txt"

     params:
         output_folder = RESULTS + "roary_st72"

     threads: 27

     shell:
         """
         rm -rf {params.output_folder}
         roary -e -n -p {threads} -f {params.output_folder} {input.gff_files}
         """

rule run_roary_ST_8:
     input:
          gff_files = [PROCESS + "prokka_STS_staphophia/" + accession + "/" + accession + ".gff" for accession in staphophia_ST_genomes["8"]] + get_gff_study(PROCESS + "prokka/", RESULTS + "MLST_results.tsv", "8")

     output:
         roary_output = RESULTS + "roary_st8/summary_statistics.txt"

     params:
         output_folder = RESULTS + "roary_st8"

     threads: 27

     shell:
         """
         rm -rf {params.output_folder}
         roary -e -n -p {threads} -f {params.output_folder} {input.gff_files}
         """

rule run_roary_ST_105:
     input:
          gff_files = [PROCESS + "prokka_STS_staphophia/" + accession + "/" + accession + ".gff" for accession in staphophia_ST_genomes["105"]] + get_gff_study(PROCESS + "prokka/", RESULTS + "MLST_results.tsv", "105")

     output:
         roary_output = RESULTS + "roary_st105/summary_statistics.txt"

     params:
         output_folder = RESULTS + "roary_st105"

     threads: 27

     shell:
         """
         rm -rf {params.output_folder}
         roary -e -n -p {threads} -f {params.output_folder} {input.gff_files}
         """

rule run_roary_ST_30:
     input:
          gff_files = [PROCESS + "prokka_STS_staphophia/" + accession + "/" + accession + ".gff" for accession in staphophia_ST_genomes["30"]] + get_gff_study(PROCESS + "prokka/", RESULTS + "MLST_results.tsv", "30")

     output:
         roary_output = RESULTS + "roary_st30/summary_statistics.txt"

     params:
         output_folder = RESULTS + "roary_st30"

     threads: 27

     shell:
         """
         rm -rf {params.output_folder}
         roary -e -n -p {threads} -f {params.output_folder} {input.gff_files}
         """

rule run_roary_ST_5:
     input:
          gff_files = [PROCESS + "prokka_STS_staphophia/" + accession + "/" + accession + ".gff" for accession in staphophia_ST_genomes["5"]] + get_gff_study(PROCESS + "prokka/", RESULTS + "MLST_results.tsv", "5") + get_gff_study(PROCESS + "prokka_references/", RESULTS + "Chile_MLST_results.tsv", "5")

     output:
         roary_output = RESULTS + "roary_st5/summary_statistics.txt"

     params:
         output_folder = RESULTS + "roary_st5"

     threads: 27

     shell:
         """
         rm -rf {params.output_folder}
         roary -e -n -p {threads} -f {params.output_folder} {input.gff_files}
         """

