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

=pod

=head1 NAME

GeneSymbolClassifier_production_conf

=head1 DESCRIPTION

eHive pipeline to assign gene symbols to the protein coding gene sequences in an Ensembl core database using a neural network gene symbol classifier.

=cut

package GeneSymbolClassifier_production_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Production::Pipeline::PipeConfig::Base_conf');

use Bio::EnsEMBL::Registry;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Hive::Version 2.5;
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

use File::Spec::Functions qw(catdir catfile);

sub default_options {
  my ($self) = @_;

  return {
    %{$self->SUPER::default_options},

    species      => [],
    taxon        => [],
    division     => [],
    run_all      => 1,
    antispecies  => [],
    antitaxons   => [],
    meta_filters => {},

    history_file => undef,

    pipeline_name => 'GeneSymbolClassifier_production',

    loading_threshold => 0.7,

    old_server_uri => '',
    ensembl_root_dir => $ENV{ENSEMBL_ROOT_DIR},
    scripts_directory => catdir($self->o('ensembl_root_dir'), 'ensembl-genes', 'pipelines', 'gene_symbol_classifier'),
    gsc_image_version => '0.12.1',
    filter_gsc_image_version => '0.3.0',
    classifier_filename => 'mlp_10_min_frequency_2022-01-29_03.08.32.ckpt',
    classifier_root_dir => '/hps/software/users/ensembl/genebuild/gene_symbol_classifier',
    classifier_img_dir => catdir($self->o('classifier_root_dir'), 'singularity'),
    classifier_data_dir => catdir($self->o('classifier_root_dir'), 'data'),
    xref_primary_accession_file => catfile($self->o('classifier_data_dir'), 'display_name_dbprimary_acc_105.dat'),
    # These images can be replace with their docker equivalent 'docker://ensemblorg/gene_symbol_classifier:0.8.2'
    gsc_image => catfile($self->o('classifier_img_dir'), 'gene_symbol_classifier_'.$self->o('gsc_image_version').'.sif'),
    filter_gsc_image => catfile($self->o('classifier_img_dir'), 'gene_symbol_classifier_filter_'.$self->o('filter_gsc_image_version').'.sif'),
  };
}

# Implicit parameter propagation throughout the pipeline.
sub hive_meta_table {
  my ($self) = @_;

  return {
    %{$self->SUPER::hive_meta_table},
    'hive_use_param_stack'  => 1,
  };
}

sub pipeline_wide_parameters {
 my ($self) = @_;

 return {
   %{$self->SUPER::pipeline_wide_parameters},
   pipeline_dir => $self->o('pipeline_dir'),
 };
}

sub pipeline_create_commands {
    my ($self) = @_;

    my $gsc_working_dir = $self->o('pipeline_dir');

    return [
        @{ $self->SUPER::pipeline_create_commands },

        "mkdir --parents --verbose $gsc_working_dir",
    ];
}

sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name      => 'InitializePipeline',
      -module          => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids       => [ {} ],
      -max_retry_count => 1,
      -flow_into       => {
        '1->A' => ['get_core_db_names'],
        'A->1' => ['Notify'],
      },
      -rc_name         => 'default',
    },

    {
      -logic_name      => 'get_core_db_names',
      -module          => 'Bio::EnsEMBL::Production::Pipeline::Common::SpeciesFactory',
      -parameters      => {
        species      => $self->o('species'),
        taxon     => $self->o('taxon'),
        division     => $self->o('division'),
        run_all      => $self->o('run_all'),
        antispecies  => $self->o('antispecies'),
        antitaxons  => $self->o('antitaxons'),
        meta_filters => $self->o('meta_filters'),
        run_all => 1,
      },
      -max_retry_count => 1,
      -flow_into       => {
        '2' => ['dump_protein_sequences'],
      },
      -rc_name         => 'default',
    },

    {
      # input: Ensembl core db
      # output: FASTA file with protein coding genes canonical sequences
      -logic_name => 'dump_protein_sequences',
      -comment    => 'Retrieve the protein coding gene sequences from a Ensembl core database and store them as a FASTA file. The analysis is auto-seeded with a job for the target core database.',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('scripts_directory').'/dump_protein_sequences.pl --group #group# --species #species# --registry #registry# --output_file #pipeline_dir#/#species#_protein_sequences.fa',
        registry => $self->o('registry'),
      },
      -max_retry_count => 1,
      -hive_capacity   => 50,
      -batch_size      => 1,
      -rc_name    => 'default',
      -flow_into  => {
        1 => 'assign_gene_symbols',
      },
    },

    {
      # input: FASTA file with protein coding genes canonical sequences
      # output: CSV file with symbol assignments and metadata
      -logic_name => 'assign_gene_symbols',
      -comment    => 'Use a gene symbol classifier neural network to assign gene symbols to protein sequences in the FASTA file and save the assignments to a CSV file.',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'singularity run --bind '.$self->o('classifier_data_dir').':/app/checkpoints --bind #pipeline_dir#:/app/data '.$self->o('gsc_image').' --checkpoint /app/checkpoints/'.$self->o('classifier_filename').q( --sequences_fasta /app/data/#species#_protein_sequences.fa --scientific_name "#expr(ucfirst(join(' ', grep {$_ !~ /gca\d+/} split('_', #species#))))expr#"),
      },
      -max_retry_count => 1,
      -hive_capacity   => 50,
      -batch_size      => 60,
      -rc_name    => 'normal',
      -flow_into  => {
        1 => 'filter_assignments',
      },
    },

    {
      # input: CSV file with symbol assignments and metadata
      # output: CSV file with assignments to be loaded to the Ensembl core db
      -logic_name => 'filter_assignments',
      -comment    => 'Filter assignments using the threshold probability and save them to a separate CSV file.',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'singularity run --bind #pipeline_dir#:/app/data '.$self->o('filter_gsc_image').' --symbol_assignments /app/data/#species#_protein_sequences_symbols.csv --threshold '.$self->o('loading_threshold'),
      },
      -max_retry_count => 1,
      -hive_capacity   => 50,
      -batch_size      => 100,
      -rc_name    => 'normal',
      -flow_into  => {
        1 => 'load_gene_symbols',
      },
    },

    {
      # input: CSV file with assignments to be loaded to the Ensembl core db
      # output: gene symbols loaded to the Ensembl core db
      -logic_name => 'load_gene_symbols',
      -comment    => 'Read gene symbols assignments from a CSV file and load them to the Ensembl core database.',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('scripts_directory').'/load_gene_symbols.pl --species #species# --group #group# --registry #registry# --symbol_assignments #pipeline_dir#/#species#_protein_sequences_symbols_filtered.csv --primary_ids_file '.$self->o('xref_primary_accession_file').' --program_version '.$self->o('gsc_image_version'),
        registry => $self->o('registry'),
      },
      -max_retry_count => 1,
      -hive_capacity   => 50,
      -batch_size      => 30,
      -rc_name    => 'default',
      -flow_into  => {
        1 => 'RunDataChecks',
      },
    },

    {
      -logic_name      => 'RunDataChecks',
      -module          => 'Bio::EnsEMBL::DataCheck::Pipeline::RunDataChecks',
      -parameters      => {
        datacheck_groups => ['xref'],
        history_file    => $self->o('history_file'),
        failures_fatal  => 1,
        output_file     => catfile('#pipeline_dir#', '#species#_dc.log'),
        registry_file   => $self->o('registry'),
        old_server_uri  => $self->o('old_server_uri'),
      },
      -max_retry_count => 1,
      -hive_capacity   => 50,
      -batch_size      => 10,
      -rc_name         => 'default',
    },

    {
      -logic_name => 'Notify',
      -module     => 'Bio::EnsEMBL::Production::Pipeline::Production::EmailSummaryCore',
      -parameters => {
        email   => $self->o('email'),
        subject => $self->o('pipeline_name').' has finished',
      },
      -rc_name    => 'default',
    },

  ];
}

sub resource_classes {
  my ($self) = @_;

  return {
    %{$self->SUPER::resource_classes},
    'normal' => { 'LSF' => ['-q production', $self->beekeeper_extra_cmdline_options]},
    'default'  => { 'LSF' => ['-q production -M 500 -R "rusage[mem=500]"', $self->beekeeper_extra_cmdline_options]},
  }
}

1;
