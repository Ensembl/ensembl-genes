#!/usr/bin/env perl


=pod

=head1 SYNOPSIS

# load gene symbols:

> perl load_gene_symbols.pl --db_host <database server hostname> --db_port <database server port> --db_name <database name> --username <database user> --password <database user password> --symbol_assignments <gene symbol assignments file path>

=head1 DESCRIPTION

Load gene symbol assignments from a file to a Ensembl core database.

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

use Bio::EnsEMBL::DBSQL::DBAdaptor;


my $db_host;
my $db_port;

my $username;
my $password;

my $db_name;
my $symbol_assignments;


# parse command line arguments
my $options = GetOptions (
    "db_host=s" => \$db_host,
    "db_port=s" => \$db_port,
    "username=s" => \$username,
    "password=s" => \$password,
    "db_name=s" => \$db_name,
    "username=s" => \$username,
    "password=s" => \$password,
    "symbol_assignments=s" => \$symbol_assignments
);


# create a database adaptor
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $db_host,
    -port   => $db_port,
    -user   => $username,
    -pass   => $password,
    -dbname => $db_name
);


# generate a gene adaptor for the database
my $gene_adaptor = $db->get_GeneAdaptor();


# open the gene symbol assignments TSV file
open(my $symbol_assignments_file, "<", $symbol_assignments);

# read gene symbol assignments from the TSV file and load them to the core database
while (my $line = <$symbol_assignments_file>) {
    print("\n");

    # remove newline (line feed and carriage return) from the line
    $line =~ s/\n//g;
    $line =~ s/\r//g;

    # get TSV field values
    my @fields = split(/\t/, $line);
    #dump(@fields);
    my $stable_id = $fields[0];
    my $symbol = $fields[1];
    my $probability = $fields[2];

    # skip TSV header line
    if ($fields[0] eq "stable_id") {
        next;
    }

    # generate gene object for gene stable_id
    my $gene = $gene_adaptor->fetch_by_stable_id($stable_id);

    # skip gene if display_xref is already defined
    if ($gene->display_xref()) {
      next;
    }

    dump($stable_id);
    dump($symbol);
    dump($probability);

    my $display_xref = $gene->display_xref();
    dump($display_xref);

    last;
}

close $symbol_assignments_file;
