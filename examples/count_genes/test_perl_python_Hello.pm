package test_perl_python_Hello;

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

'output_path'                  => '/hps/nobackup2/production/ensembl/kbillis/production/tasks/fix_genes_with_introns/fergals_pipeline/',
'production_name'              => 'sarcophilus_harrisii', # The production name, including any modifiers
'species_name'                 => 'sarcophilus_harrisii', # e.g. mus_musculus
'taxon_id'                     => '9305',
'genus_taxon_id'               => '' || $self->o('taxon_id'),


'farm_user_name'        => 'kbillis', # for output_db prefix
'enscode_root_dir'      => '/nfs/production/panda/ensembl/kbillis/enscode_2019_11/enscode/', # git repo checkouts
'user_r'                => 'ensro',
'user'                  => 'ensadmin',
'password'              => 'ensembl',

'dna_dbname'            => 'kbillis_sarcophilus_harrisii_core_100_tmp',
'dna_db_server'         => $ENV{GBS1},
'dna_db_port'           => $ENV{GBP1},


'pipeline_name'         => 'kbillis_python_perl',
'pipe_db_server'        => $ENV{GBS4},
'port'                  => $ENV{GBP4},


'dna_dbname_DELETE'            => 'kbillis_sarcophilus_harrisii_core_100_tmp_DELETEME',



'base_blast_db_path'     => $ENV{BLASTDB_DIR},
'uniprot_version'        => 'uniprot_2018_07',
'protein_blast_db'       => '' || catfile($self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'PE12_vertebrata'),
'protein_blast_index'    => '' || catdir($self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'PE12_vertebrata_index'),

'blast_type' => 'ncbi', # It can be 'ncbi', 'wu', or 'legacy_ncbi'
'uniprot_blast_exe_path' => catfile($self->o('binary_base'), 'blastp'),
'use_threads'           => 3,



'dna_db' => {
              -dbname => $self->o('dna_dbname'),
              -host   => $self->o('dna_db_server'),
              -port   => $self->o('dna_db_port'),
              -user   => $self->o('user_r'),
              -driver => $self->o('hive_driver'),
            },





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
      	  # before you need:  
      	  # export PYTHONPATH="/nfs/production/panda/ensembl/kbillis/python_in_hive/ensembl-genes/examples/count_genes/:/nfs/production/panda/ensembl/kbillis/code/python/modules/:/nfs/production/panda/ensembl/kbillis/code/python/modules/lib/python2.7/site-packages/"
      	  -logic_name => 'take_python_part',
          -module     => 'examples.count_genes.hello_world',
          -language   => 'python3',
          # -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
          -analysis_capacity  =>  2,  # use per-analysis limiter
          -wait_for => ['Check_perl_is'],
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
