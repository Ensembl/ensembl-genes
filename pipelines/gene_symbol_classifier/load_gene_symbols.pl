#!/usr/bin/env perl


=pod

=head1 SYNOPSIS

# load gene symbols:

> perl load_gene_symbols.pl --db_host <database server hostname> --db_port <database server port> --db_name <database name> --username <database user> --password <database user password> --symbol_assignments <gene symbol assignments file path>

=head1 DESCRIPTION

Load gene symbol assignments from a CSV file to a Ensembl core database.

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


use strict;
use warnings;
use feature "say";

use Getopt::Long;

use Data::Dump qw(dump);

use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


my $db_host;
my $db_port;

my $username;
my $password;

my $db_name;
my $symbol_assignments;


my %db_display_name_to_db_name = (
    "Clone-based (Ensembl) gene" => "Clone_based_ensembl_gene",
    "NCBI gene (formerly Entrezgene)" => "EntrezGene",
    "FlyBase gene name" => "FlyBaseName_gene",
    "FlyBase annotation" => "flybase_annotation_id",
    "HGNC Symbol" => "HGNC",
    "MGI Symbol" => "MGI",
    "RGD Symbol" => "RGD",
    "SGD gene name" => "SGD_GENE",
    "UniProtKB Gene Name" => "Uniprot_gn",
    "VGNC Symbol" => "VGNC",
    "WormBase Gene Sequence-name" => "wormbase_gseqname",
    "WormBase Locus" => "wormbase_locus",
    "Xenbase" => "Xenbase",
    "ZFIN" => "ZFIN_ID",
);


# parse command line arguments
my $options = GetOptions (
    "db_host=s" => \$db_host,
    "db_port=s" => \$db_port,
    "username=s" => \$username,
    "password=s" => \$password,
    "db_name=s" => \$db_name,
    "username=s" => \$username,
    "password=s" => \$password,
    "symbol_assignments=s" => \$symbol_assignments,
);


# create a database adaptor
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $db_host,
    -port   => $db_port,
    -user   => $username,
    -pass   => $password,
    -dbname => $db_name,
);


# generate a gene adaptor
my $gene_adaptor = $db->get_GeneAdaptor();

# generate a DBEntry adaptor for Xref
my $xref_adaptor = $db->get_DBEntryAdaptor();

# open the gene symbol assignments CSV file
open(my $symbol_assignments_file, "<", $symbol_assignments) or die "Could not open file $symbol_assignments: $!";

# read gene symbol assignments from the CSV file and load them to the core database
while (my $line = <$symbol_assignments_file>) {
    # remove newline (line feed and carriage return) from the line
    $line =~ s/\n//g;
    $line =~ s/\r//g;

    # get CSV field values
    my @fields = split(/\t/, $line);
    my $stable_id = $fields[0];
    my $symbol = $fields[1];
    my $probability = $fields[2];
    my $symbol_description = $fields[3];
    my $symbol_source = $fields[4];

    # skip CSV header line
    if ($fields[0] eq "stable_id") {
        next;
    }

    # generate gene object for gene stable_id
    my $gene = $gene_adaptor->fetch_by_stable_id($stable_id);

    # skip gene if display_xref is already set
    my $display_xref = $gene->display_xref();
    if ($display_xref) {
        say("gene ".$stable_id." already has an Xref, symbol: ".$display_xref->display_id().", skipping");
        next;
    }

    # generate Xref object
    my $score = 100 * sprintf("%.4f", $probability);
    my $assignment_description = "Ensembl NN prediction with score ".$score."%";
    my $xref_dbname = $db_display_name_to_db_name{$symbol_source};
    my $xref_object = Bio::EnsEMBL::DBEntry->new(
        -primary_id => $symbol,
        -display_id => $symbol,
        -dbname => $xref_dbname,
        -description => $symbol_description,
    );

    my $gene_description = $symbol_description." [".$assignment_description."]";

    # store the Xref to the database and add it to the gene
    $xref_adaptor->store($xref_object, $gene->dbID, "Gene", 1);
    $gene->display_xref($xref_object);
    $gene->description($gene_description);
    $gene_adaptor->update($gene);

    say("added symbol ".$symbol.", ".$gene_description." for gene ".$stable_id);
}

close $symbol_assignments_file;
