=pod

=head1 NAME

gene_symbol_classifier_conf;

=head1 SYNOPSIS

# pipeline initialization:

> init_pipeline.pl gene_symbol_classifier_conf --pipe_db_server <pipeline MySQL server hostname> --pipe_db_port <pipeline MySQL server port> --user <username> --password <password> --user_r ensro --pipe_db_name <pipeline database name> --core_db_server_host <core database server hostname> --core_db_server_port <core database server port> --core_db_name <core database name> --annotation_data_directory <annotation data directory> --singularity_image <singularity image path> --classifier_directory <classifier checkpoint directory> --classifier_filename <classifier checkpoint filename> --scientific_name <assembly scientific name>

=head1 DESCRIPTION

eHive pipeline to assign gene symbols to the protein coding gene sequences in an Ensembl core database using a neural network gene symbol classifier.

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


package gene_symbol_classifier_conf;


use strict;
use warnings;

use parent ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');


sub default_options {
    my ($self) = @_;

    # TODO
    # replace when incorporated to the main Genebuild annotation pipeline
    #my $annotation_data_directory = $self->o('annotation_data_directory');
    #my $gsc_data_directory = "${annotation_data_directory}/gene_symbol_classifier";
    my $gsc_data_directory = $self->o('annotation_data_directory');

    my $core_db_name = $self->o('core_db_name');

    my $protein_sequences_fasta_path = "${gsc_data_directory}/${core_db_name}_protein_sequences.fa";
    my $gene_symbols_csv_path = "${gsc_data_directory}/${core_db_name}_protein_sequences_symbols.csv";

    return {
        # inherit from the base class
        %{ $self->SUPER::default_options() },

        'gsc_data_directory' => $gsc_data_directory,
        'protein_sequences_fasta_path' => $protein_sequences_fasta_path,
        'gene_symbols_csv_path' => $gene_symbols_csv_path,

        'pipeline_db' => {
            -driver => 'mysql',
            -host   => $self->o('pipe_db_server'),
            -port   => $self->o('pipe_db_port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -dbname => $self->o('pipe_db_name'),
        },
    };
}


sub pipeline_create_commands {
    my ($self) = @_;

    my $gsc_data_directory = $self->o('gsc_data_directory');

    return [
        @{ $self->SUPER::pipeline_create_commands },

        "mkdir --parents --verbose $gsc_data_directory",
    ];
}


sub pipeline_analyses {
    my ($self) = @_;

    return [
        {
            # input: Ensembl core db
            # output: FASTA file
            -logic_name => 'dump_protein_sequences',
            -comment    => 'Retrieve the protein coding gene sequences from a Ensembl core database and store them as a FASTA file. The analysis is auto-seeded with a job for the target core database.',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -input_ids  => [ {} ],
            -parameters => {
                'cmd' => 'if [ ! -e "'.$self->o('protein_sequences_fasta_path').'" ]; then perl dump_protein_sequences.pl --db_host '.$self->o('core_db_server_host').' --db_port '.$self->o('core_db_server_port').' --db_name '.$self->o('core_db_name').' --output_file '.$self->o('protein_sequences_fasta_path').'; fi',
            },
            -rc_name    => 'default',
            -flow_into  => {
                1 => 'assign_gene_symbols',
            },
        },

        {
            # input: FASTA file
            # output: CSV file
            -logic_name => 'assign_gene_symbols',
            -comment    => 'Use a gene symbol classifier neural network to assign gene symbols to protein sequences in the FASTA file and save the assignments to a CSV file.',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd' => 'if [ ! -e "'.$self->o('gene_symbols_csv_path').'" ]; then singularity run --bind '.$self->o('classifier_directory').':/app/checkpoints --bind '.$self->o('gsc_data_directory').':/app/data '.$self->o('singularity_image').' --checkpoint /app/checkpoints/'.$self->o('classifier_filename').' --sequences_fasta /app/data/'.$self->o('core_db_name').'_protein_sequences.fa --scientific_name "'.$self->o('scientific_name').'"; fi',
            },
            -rc_name    => 'default',
            -flow_into  => {
                1 => 'load_gene_symbols',
            },
        },

        {
            # input: CSV file
            # output: gene symbols added to Ensembl core db
            -logic_name => 'load_gene_symbols',
            -comment    => 'Read gene symbols assignments from a CSV file and load them to the Ensembl core database.',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd' => 'perl load_gene_symbols.pl --db_host '.$self->o('core_db_server_host').' --db_port '.$self->o('core_db_server_port').' --db_name '.$self->o('core_db_name').' --username '.$self->o('user').' --password '.$self->o('password').' --symbol_assignments '.$self->o('gene_symbols_csv_path'),
            },
            -rc_name    => 'default',
        },
    ];
}


sub resource_classes {
    my $self = shift;

    return {
        'default' => {
            'LSF' => $self->lsf_resource_builder(
                'production',
                900,
                [
                    $self->default_options->{'pipe_db_server'},
                    $self->default_options->{'dna_db_server'}
                ],
                [ $self->default_options->{'num_tokens'} ]
            )
        },
    };
}


1;
