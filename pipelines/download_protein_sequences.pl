#!/usr/bin/env perl


=pod

=head1 SYNOPSIS

# download protein coding gene sequences:

> perl download_protein_sequences.pl --db_host <database server hostname> --db_port <database server port> --db_name <database name> --output_file <output file path>

=head1 DESCRIPTION

Retrieve the protein coding gene sequences from a Ensembl core database and store them as a FASTA file.

=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


use warnings;
use strict;
use feature "say";

use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;


my $db_host;
my $db_port;

my $username = "ensro";
my $password;

my $db_name;
my $output_file;


# parse command line arguments
my $options = GetOptions (
    "db_host=s" => \$db_host,
    "db_port=s" => \$db_port,
    "username=s" => \$username,
    "password=s" => \$password,
    "db_name=s" => \$db_name,
    "output_file=s" => \$output_file
);


# create a database adaptor
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $db_host,
    -port   => $db_port,
    -user   => $username,
    -pass   => $password,
    -dbname => $db_name
);


my $slice_adaptor = $db->get_SliceAdaptor();
my $meta_adaptor = $db->get_MetaContainer();

my $coordinate_system = "toplevel";
my $slices = $slice_adaptor->fetch_all($coordinate_system);
my $production_name = $meta_adaptor->get_production_name;
say $production_name." genome assembly";

# open output FASTA file for writing
open(my $sequences_fasta_file, ">", $output_file);

# iterate slices and save protein coding genes stable ids and sequences to the output file
foreach my $slice (@$slices) {
    # get genes that overlap this slice
    my $genes = $slice->get_all_Genes();

    say "processing seq_region ".$slice->seq_region_name." with ".scalar(@$genes)." total genes";

    foreach my $gene (@$genes) {
        # skip non-protein coding genes
        my $biotype = $gene->biotype;
        unless ($biotype eq "protein_coding") {
            next;
        }

        # get the canonical transcript of the gene
        my $transcript = $gene->canonical_transcript;
        unless ($transcript) {
            die "no canonical transcript for gene".$gene->stable_id;
        }

        my $gene_stable_id = $gene->stable_id;
        my $transcript_stable_id = $transcript->stable_id;

        my $translation = $transcript->translation;
        my $translation_stable_id = $translation->stable_id;
        my $protein_sequence = $translation->seq;

        say $sequences_fasta_file ">$gene_stable_id translation_stable_id:$translation_stable_id transcript_stable_id:$transcript_stable_id";
        say $sequences_fasta_file $protein_sequence;
    }
}

close $sequences_fasta_file;
