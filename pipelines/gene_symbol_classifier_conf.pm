=pod

=head1 NAME

    gene_symbol_classifier_conf;

=head1 SYNOPSIS

    # pipeline initialization:

    > init_pipeline.pl gene_symbol_classifier_conf -pipe_db_server <pipeline MySQL server host> -pipe_db_port <pipeline MySQL server port> -user <username> -password <password> -user_r ensro -core_db <core db name> -core_db_server_host <core db server host> -core_db_server_port <core db server port> -annotation_data_directory <annotation data directory>

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
    #my $gene_symbol_classifier_directory = "${annotation_data_directory}/gene_symbol_classifier";
    my $gene_symbol_classifier_directory = $self->o('annotation_data_directory');

    my $protein_sequences_fasta_path = "${gene_symbol_classifier_directory}/gene_sequences.fasta";
    my $gene_symbols_tsv_path = "${gene_symbol_classifier_directory}/gene_symbols.tsv";

    return {
        # inherit from the base class
        %{ $self->SUPER::default_options() },

        'protein_sequences_fasta_path' => $protein_sequences_fasta_path,
        'gene_symbols_tsv_path' => $gene_symbols_tsv_path,

        'pipeline_db' => {
            -driver => 'mysql',
            -host   => $self->o('pipe_db_server'),
            -port   => $self->o('pipe_db_port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -dbname => 'gsc_testing',
        },
    };
}


sub pipeline_create_commands {
    my ($self) = @_;

    # TODO
    # replace when incorporated to the main Genebuild annotation pipeline
    #my $annotation_data_directory = $self->o('annotation_data_directory');
    #my $gene_symbol_classifier_directory = "${annotation_data_directory}/gene_symbol_classifier";
    my $gene_symbol_classifier_directory = $self->o('annotation_data_directory');

    return [
        @{ $self->SUPER::pipeline_create_commands },

        "mkdir --parents --verbose $gene_symbol_classifier_directory",
    ];
}


sub pipeline_analyses {
    my ($self) = @_;

    return [
        {
            # input: Ensembl core db
            # output: FASTA file
            -logic_name => 'dump_fasta',
            -comment    => 'Retrieve the protein coding gene sequences from a Ensembl core database and store them as a FASTA file. The analysis is auto-seeded with a job for the target core database.',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -input_ids  => [
                {
                    'core_db' => $self->o('core_db'),
                    'core_db_server_host' => $self->o('core_db_server_host'),
                    'core_db_server_port' => $self->o('core_db_server_port'),
                }
            ],
            -parameters => {
                'cmd' => 'echo "dump_fasta analysis"',
            },
            -rc_name    => 'default',
            -flow_into  => {
                1 => 'run_classifier',
            },
        },

        {
            # input: FASTA file
            # output: TSV file
            -logic_name => 'run_classifier',
            -comment    => 'Use a gene symbol classifier neural network to assign gene symbols to protein sequences in the FASTA file and save the assignments to a TSV file.',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd' => 'echo "run_classifier analysis"',
            },
            -rc_name    => 'default',
            -flow_into  => {
                1 => 'load_gene_symbols',
            },
        },

        {
            # input: TSV file
            # output: gene symbols added to Ensembl core db
            -logic_name => 'load_gene_symbols',
            -comment    => 'Read gene symbols assignments from a TSV file and load them to the Ensembl core database.',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd' => 'echo "load_gene_symbols analysis"',
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
                'production-rh74',
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
