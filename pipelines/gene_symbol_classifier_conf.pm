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

    return {
        # inherit from the base class
        %{ $self->SUPER::default_options() },

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

    return [ @{ $self->SUPER::pipeline_create_commands }, ];
}


sub pipeline_analyses {
    my ($self) = @_;

    return [
        {
            # input: genome assembly core db
            # output: FASTA file
            -logic_name => 'dump_fasta',
            -comment    => 'retrieve the protein coding gene sequences from a genome assembly core database and store them as a FASTA file',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd' => 'echo "dump_fasta analysis"',
            },
            -rc_name   => 'default',
            -input_ids => [ {} ],
            -flow_into => { '1' => 'run_classifier' }
        },

        {
            # input: FASTA file
            # output: TSV file
            -logic_name => 'run_classifier',
            -comment    => 'use a gene symbol classifier neural network to assign symbols to gene protein sequences in a FASTA file and save them in a TSV file',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd' => 'echo "run_classifier analysis"',
            },
            -rc_name   => 'default',
            -flow_into => { '1' => 'load_gene_symbols' }
        },

        {
            # input: TSV file
            # output: genome assembly core db
            -logic_name => 'load_gene_symbols',
            -comment    => 'read a TSV file with gene symbols and load them to the genome assembly core database',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd' => 'echo "load_gene_symbols analysis"',
            },
            -rc_name => 'default',
        },
    ];
}


sub resource_classes {
    my $self = shift;

    return {
        'default' => {
            LSF => $self->lsf_resource_builder(
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
