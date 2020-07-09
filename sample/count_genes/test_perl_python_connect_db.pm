package test_perl_python_connect_db;

use strict;
use warnings;
use File::Spec::Functions;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');




sub default_options {
  my ($self) = @_;
  return {
    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },

   'pipeline_name'         => 'kbillis_python_perl',
   'pipe_db_server'        => $ENV{GBS4},
   'port'                  => $ENV{GBP4},

    'user_r'                => 'ensro',
    'user'                  => 'ensadmin',
    'password'              => 'ensembl',

    'dna_dbname'            => 'kbillis_sarcophilus_harrisii_core_100_tmp',
    'dna_db_server'         => $ENV{GBS1},
    'dna_db_port'           => $ENV{GBP1},


  
   'use_threads'           => 3,

 } # end return
} # end default_options


sub pipeline_create_commands {
    my ($self) = @_;
    return [
    # inheriting database and hive tables' creation
	    @{$self->SUPER::pipeline_create_commands},
    ];
} # end pipeline_create_commands


sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
	  %{$self->SUPER::pipeline_wide_parameters},
    
        'take_time'     => 1,
    
  }
}

sub hive_meta_table {
    my ($self) = @_;
    return {
        %{$self->SUPER::hive_meta_table},       # here we inherit anything from the base class

        'hive_use_param_stack'  => 1,           # switch on the new param_stack mechanism
    };
}

sub pipeline_analyses {
  my ($self) = @_;
  return [
  

      {
        -logic_name => 'Check_perl_is',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ `perl -version | wc -c` -ne "0" ]; then echo "Hey Perl"; else echo "PERL?" ;fi',
                         # return_codes_2_branches => {'42' => 2},
                       },
        -flow_into  => {
          1 => ['take_python_part'],
        },
        -rc_name => 'default',
        -input_ids => [{
                 'hello' => 'Yo!'         	
        }],
      },

      {   
      	  -logic_name => 'take_python_part',
          -module     => 'sample.count_genes.hello_world',
          -language   => 'python3',
          # -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
          -analysis_capacity  =>  2,  # use per-analysis limiter
          -rc_name => 'default',
          -flow_into  => {
            1 => ['python_connect_db'],
          },
          
      },

      {   
      	  -logic_name => 'python_connect_db',
          -module     => 'sample.count_genes.connect_to_mysql_db',
          -language   => 'python3',
          -input_ids => [
                { 'db_name' => 'kbillis_sarcophilus_harrisii_core_100', 'port' => '4532' } 
                ], 
      
          # -db_connection_list => '["mysql-ens-genebuild-prod-6","4532", "ensro", "kbillis_sarcophilus_harrisii_core_100"]'

          # -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
          -wait_for => ['take_python_part'],
          -rc_name => 'default',
          -analysis_capacity  =>  2,  # use per-analysis limiter
      },

  ]
} # end pipeline analyses


sub resource_classes {
  my $self = shift;

  return {
    'default' => { 'LSF' => $self->lsf_resource_builder('production-rh74',500) },
    '1GB' =>  { 'LSF' => $self->lsf_resource_builder('production-rh74',1000) },
    '2GB' =>  { 'LSF' => $self->lsf_resource_builder('production-rh74',2000) },
    '5GB' =>  { 'LSF' => $self->lsf_resource_builder('production-rh74',5000) },
    '10GB' => { 'LSF' => $self->lsf_resource_builder('production-rh74',10000) },
    '15GB' => { 'LSF' => $self->lsf_resource_builder('production-rh74',15000) },
    '20GB' => { 'LSF' => $self->lsf_resource_builder('production-rh74',20000) },
    '25GB' => { 'LSF' => $self->lsf_resource_builder('production-rh74',25000) },
    '30GB' => { 'LSF' => $self->lsf_resource_builder('production-rh74',30000) },
    'blast' => { LSF => $self->lsf_resource_builder('production-rh74', 2000, undef, undef, ($self->default_options->{'use_threads'}+1))},
    'blast10GB' => { LSF => $self->lsf_resource_builder('production-rh74', 10000, undef, undef, ($self->default_options->{'use_threads'}+1))},
  }
} # end resource_classes

1;
