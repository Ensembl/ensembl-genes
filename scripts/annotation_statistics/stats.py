import gffutils, argparse, sys

parser = argparse.ArgumentParser(description="Calculate stats from a GFF3 file")
parser.add_argument(
    "-f",
    "--gff3_file",
    help="Name of GFF3 file",
    required=True,
)
args = parser.parse_args()
gff3_file = args.gff3_file

file_prefix = "-".join(gff3_file.split("-", 2)[:2])

# Create a database from the GFF3 file
db = gffutils.create_db(
    gff3_file, ":memory:", id_spec="ID", merge_strategy="create_unique"
)

# Initialize counters
num_coding_genes = 0
num_noncoding_genes = 0
num_transcripts_per_coding_gene = 0
num_exons_per_transcript = 0
coding_seq_length = 0
coding_exon_length = 0
coding_intron_length = 0

# Iterate through the features in the database
for gene in db.features_of_type("gene"):
    if "biotype" in gene.attributes:
        biotype = gene.attributes["biotype"][0]
        if biotype == "protein_coding":
            num_coding_genes += 1
            transcripts = list(db.children(gene, featuretype="mRNA"))
            num_transcripts_per_coding_gene += len(transcripts)
            for transcript in transcripts:
                exons = list(db.children(transcript, featuretype="exon"))
                num_exons_per_transcript += len(exons)
                for exon in exons:
                    coding_exon_length += exon.end - exon.start
                # come back to intron code - feature "intron" does not exist, will need to calculate as length between exons
                #                if exons:
                #                    introns = list(db.children(transcript, featuretype='intron'))
                #                    for intron in introns:
                #                        coding_intron_length += intron.end - intron.start
                cdses = list(db.children(transcript, featuretype="CDS"))
                for cds in cdses:
                    coding_seq_length += cds.end - cds.start

for gene in db.features_of_type("ncRNA_gene"):
    num_noncoding_genes += 1

# Calculate averages
avg_transcripts_per_coding_gene = num_transcripts_per_coding_gene / num_coding_genes
avg_exons_per_transcript = num_exons_per_transcript / num_transcripts_per_coding_gene
avg_coding_seq_length = coding_seq_length / num_transcripts_per_coding_gene
avg_coding_exon_length = coding_exon_length / num_exons_per_transcript
# avg_coding_intron_length = coding_intron_length / (num_exons_per_transcript - num_transcripts_per_coding_gene)

# Print the results
output_file = open(file_prefix + "_annotation_statistics.txt", "w")
output = (
    "PCG_NUM\t"
    + str(num_coding_genes)
    + "\nNCG_NUM\t"
    + str(num_noncoding_genes)
    + "\nAVG_TR_PER_GENE\t"
    + str(avg_transcripts_per_coding_gene)
    + "\nAVG_EXON_PER_TR\t"
    + str(avg_exons_per_transcript)
    + "\nAVG_CDS_LEN\t"
    + str(avg_coding_seq_length)
    + "\nAVG_EXON_LEN\t"
    + str(avg_coding_exon_length)
)
# print('Average length of coding intron size:', avg_coding_intron_length)

print(output, file=output_file)
