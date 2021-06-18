
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

package Bio::EnsEMBL::Analysis::Hive::Config::gene_symbol_classifier;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

sub default_options {
  my ($self) = @_;
  return {
    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },
    'pipeline_db' => {
      -dbname => 'ws_gsc',
      -host   => 'mysql-ens-genebuild-prod-2',
      -port   => 4528,
      -user   => 'ensadmin',
      -pass   => 'ensembl',
      -driver => 'mysql',
    },
  };
}

sub pipeline_create_commands {
  my ($self) = @_;
  return [ @{ $self->SUPER::pipeline_create_commands }, ];
}

sub pipeline_analyses {
  my ($self) = @_;
  return [];
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
