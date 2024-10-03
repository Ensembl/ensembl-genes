use warnings;
use strict;
use feature 'say';
use Data::Dumper;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);
use List::MoreUtils qw(first_index);

my $coord_system = 'toplevel';
my $user = 'ensro';
my $dbname = '';
my $host = '';
my $port = '';
my $pass;
my $production_name = '';
my $assembly_names_input_file = '';
my $output_dir = '';
GetOptions(
    'host=s'   => \$host,
    'port=s'   => \$port,
    'dbname=s' => \$dbname,
    'production_name=s' => \$production_name,
    'assembly_names=s' => \$assembly_names_input_file,
    'output_dir=s' => \$output_dir
);

my @coding_headers = ('Scientific name',
                      'Coding genes',
                      'Average genomic span',
                      'Average sequence length',
                      'Average CDS length',
                      'Shortest gene',
                      'Longest gene',
                      'Total transcripts',
                      'Coding transcripts',
                      'Transcripts per gene',
                      'Coding transcripts per gene',
                      'Total exons',
                      'Total coding exons',
                      'Average exon length',
                      'Average coding exon length',
                      'Average exons per transcript',
                      'Average coding exons per coding transcript',
                      'Total introns',
                      'Average intron length',
                      );

my @non_coding_headers = ('Scientific name',
                          'Non-coding genes',
                          'Small non-coding genes',
                          'Long non-coding genes',
                          'Misc non-coding genes',
                          'Average genomic span',
                          'Average sequence length',
                          'Shortest gene',
                          'Longest gene',
                          'Total transcripts',
                          'Transcripts per gene',
                          'Total exons',
                          'Average exon length',
                          'Average exons per transcript',
                          'Total introns',
                          'Average intron length',
                          );

my @pseudogene_headers = ('Scientific name',
                          'Pseudogenes',
                          'Average genomic span',
                          'Average sequence length',
                          'Shortest gene',
                          'Longest gene',
                          'Total transcripts',
                          'Transcripts per gene',
                          'Total exons',
                          'Average exon length',
                          'Average exons per transcript',
                          'Total introns',
                          'Average intron length',
                          );


my @assembly_headers = ('Scientific name',
                        'Sex',
                        'Breed/Cultivar/Strain',
                        'Taxonomy id',
                        'Assembly name',
                        'Assembly accession',
                        'Assembly date',
                        'Contig N50',
                        'Total genome length',
                        'Total coding sequence length',
                        'Total gap length',
                        'Spanned gaps',
                        'Chromosomes',
                        'Toplevel sequences',
                        'Component sequences',
                        '% GC');

my @final_coding_stats = (\@coding_headers);
my @final_non_coding_stats = (\@non_coding_headers);
my @final_pseudogene_stats = (\@pseudogene_headers);
my @final_assembly_stats = (\@assembly_headers);

my $ftp_base = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all';

  say "Processing: ".$dbname;

  # The code below is awkward for a variety of reasons. The first is because the meta table isn't consistent
  # and some values such as the GCA or assembly date, which would seem fundamental, are not always present
  # The second issue is the collection dbs, they seem to need the species id and production name to be supplied
  # when creating the adaptor. So some code was needed to lookup the species id in addition to the production name
  # Then there are some assemblies where the assembly name in the meta table is not the same as in the assembly stats
  # file in the archives, so a mapping is needed unless we decide to change the name in the meta table
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -port    => $port,
    -user    => $user,
    -host    => $host,
    -dbname  => $dbname,
    -pass    => $pass);

  my $sth = $db->dbc->prepare("SELECT species_id FROM meta WHERE meta_key='species.production_name' AND meta_value=?");
  $sth->bind_param(1,$production_name);
  say ($production_name);
  $sth->execute();
  my ($species_id) = $sth->fetchrow_array();

  $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -port       => $port,
    -user       => $user,
    -host       => $host,
    -dbname     => $dbname,
    -pass       => $pass,
    -species    => $production_name,
    -species_id => $species_id);

  say "Species id: ".$species_id;

  my $biotype_groups = { 'coding' => [],
                         'non-coding' => [],
                         'pseudogene' => [],
                       };

  my $slice_adaptor = $db->get_SliceAdaptor();
  my $slices = $slice_adaptor->fetch_all('toplevel');
  my $gene_adaptor = $db->get_GeneAdaptor();
  my $meta_adaptor = $db->get_MetaContainerAdaptor();
  my $gca = $meta_adaptor->single_value_by_key('assembly.accession');
  if ($gca =~ /^GCF_/) {
    $gca = $meta_adaptor->single_value_by_key('assembly.alt_accession', 1);
  }

  my $assembly_name = $meta_adaptor->single_value_by_key('assembly.name');
  open(my $assembly_names_file, '<', $assembly_names_input_file);
  my @lines = <$assembly_names_file>;
  close $assembly_names_file;
  my $assembly_maps;
  foreach my $line (@lines) {
    chomp $line;
    my @entries = split / /, $line;
    my $entry_db = $entries[0];
    my $entry_asm_name = $entries[1];

    $assembly_maps->{$entry_db} = $entry_asm_name;
  }

  # For some reason the assembly name is different in meta table
  if(exists($assembly_maps->{$dbname})) {
    $assembly_name = $assembly_maps->{$dbname};
  }

  my $strain = $meta_adaptor->single_value_by_key('species.strain');

  unless($strain) {
    $strain = '';
  }

  my $scientific_name = $meta_adaptor->single_value_by_key('species.scientific_name');
  my $taxon_id = $meta_adaptor->single_value_by_key('species.taxonomy_id');

  my $stats_file = $gca.'_'.$assembly_name.'_assembly_stats.txt';
  my $stats_path = $ftp_base.'/'.substr($gca, 0, 3).'/'.substr($gca, 4, 3).'/'.substr($gca, 7, 3).'/'.substr($gca, 10, 3).'/'.$gca.'_'.$assembly_name.'/'.$stats_file;

  my $stats_cmd = "wget -P ".$output_dir." ".$stats_path;
  my $result = system($stats_cmd);
  if($result) {
    die "wget command returned a non-zero exit code. Command line used:\n".$stats_cmd;
  }

  unless(-e  $output_dir."/".$stats_file) {
    die "Could not download the report file. Path used: ".$stats_path;
  }

  my $stats_keys = {
   'total-length' => '',
   'contig-N50' => '',
   'total-gap-length' => '',
   'top-level-count' => '',
   'component-count' => '',
   'molecule-count' => '',
   'spanned-gaps' => '',
   'sex' => '',
   'date' => '',
  };

  open(IN, $output_dir."/".$stats_file);
  while(<IN>) {
    my $line = $_;
    chomp($line);
    if($line =~ /^\# Sex\:\s+(.+)/) {
      $stats_keys->{'sex'} = $1;
    }

    if($line =~ /^\# Date\:\s+([\d\-]+)/) {
      $stats_keys->{'date'} = $1;
    }

    unless($line =~ /^all/) {
      next;
    }

    my @ele = split("\t",$line);
    my $key = $ele[4];
    $ele[5] =~ /\d+/;
    my $val = $&;
    unless(exists($stats_keys->{$key})) {
      next;
    }

    $stats_keys->{$key} = $val;
  }
  close IN;

  my $coding_genes = [];
  my $non_coding_genes = [];
  my $pseudogenes = [];

  my $total_sequence_length = 0;
  my $total_coding_sequence_length = 0;
  my $total_chromosome_count = 0;
  my $total_chromosome_gene_count = 0;
  my $total_chromosome_coding_gene_count = 0;
  my $total_chromosome_non_coding_gene_count = 0;
  my $total_chromosome_pseudogene_count = 0;
  my $total_gene_count = 0;
  my $total_genomic_regions_count = scalar(@$slices);
  my $total_gc_count = 0;

  foreach my $slice (@$slices) {
    my $seq = $slice->seq();
    my $gc_count = $seq =~ tr/GgCc//;
    $seq = '';
    $total_sequence_length += $slice->length();
    $total_gc_count += $gc_count;

    my $genes = $slice->get_all_Genes();
    if($slice->has_karyotype()) {
      $total_chromosome_count++;
      $total_chromosome_gene_count += scalar(@$genes);
    }

    $total_gene_count += scalar(@$genes);
    foreach my $gene (@$genes) {
      my $biotype_group = $gene->get_Biotype->biotype_group;
      if($biotype_group eq 'coding') {
        push(@$coding_genes,$gene);
        if($slice->has_karyotype()) {
          $total_chromosome_coding_gene_count++;
        }
      } elsif($biotype_group eq 'snoncoding' || $biotype_group eq 'lnoncoding' || $biotype_group eq 'mnoncoding') {
        push(@$non_coding_genes,$gene);
        if($slice->has_karyotype()) {
          $total_chromosome_non_coding_gene_count++;
        }
      } elsif($biotype_group eq 'pseudogene') {
        push(@$pseudogenes,$gene);
        if($slice->has_karyotype()) {
          $total_chromosome_pseudogene_count++;
        }
      }
    } # End foreach my $gene
  } # End foreach my $slice

  say "Coding genes: ".scalar(@$coding_genes);
  say "Non-coding genes: ".scalar(@$non_coding_genes);
  say "Pseudogenes: ".scalar(@$pseudogenes);

  my ($coding_stats,$coding_sequence_length) = process_coding_genes($coding_genes,$scientific_name);
  $total_coding_sequence_length += $coding_sequence_length;
  my $non_coding_stats = process_non_coding_genes($non_coding_genes,$scientific_name);
  my $pseudogene_stats = process_pseudogenes($pseudogenes,$scientific_name);
  my $gc_percent = sprintf("%.2f", ($total_gc_count/$total_sequence_length) * 100);

  my $assembly_stats = [$scientific_name,
                        $stats_keys->{'sex'},
                        $strain,
                        $taxon_id,
                        $assembly_name,
                        $gca,
                        $stats_keys->{'date'},
                        $stats_keys->{'contig-N50'},
                        $stats_keys->{'total-length'},
                        $total_coding_sequence_length,
                        $stats_keys->{'total-gap-length'},
                        $stats_keys->{'spanned-gaps'},
                        $stats_keys->{'molecule-count'},
                        $stats_keys->{'top-level-count'},
                        $stats_keys->{'component-count'},
                        $gc_percent,
                        ];

  push(@final_coding_stats,$coding_stats);
  push(@final_non_coding_stats,$non_coding_stats);
  push(@final_pseudogene_stats,$pseudogene_stats);
  push(@final_assembly_stats,$assembly_stats);

  #create hashes for stats meta keys
  my %coding_stats_hash = ('genebuild.stats.coding_genes' => $$coding_stats[1],
                           'genebuild.stats.average_genomic_span' => $$coding_stats[2],
                           'genebuild.stats.average_sequence_length' => $$coding_stats[3],
                           'genebuild.stats.average_cds_length' => $$coding_stats[4],
                           'genebuild.stats.shortest_gene_length' => $$coding_stats[5],
                           'genebuild.stats.longest_gene_length' => $$coding_stats[6],
                           'genebuild.stats.total_transcripts' => $$coding_stats[7],
                           'genebuild.stats.coding_transcripts' => $$coding_stats[8],
                           'genebuild.stats.transcripts_per_gene' => $$coding_stats[9],
                           'genebuild.stats.coding_transcripts_per_gene' => $$coding_stats[10],
                           'genebuild.stats.total_exons' => $$coding_stats[11],
                           'genebuild.stats.total_coding_exons' => $$coding_stats[12],
                           'genebuild.stats.average_exon_length' => $$coding_stats[13],
                           'genebuild.stats.average_coding_exon_length' => $$coding_stats[14],
                           'genebuild.stats.average_coding_exons_per_transcript' => $$coding_stats[15],
                           'genebuild.stats.average_coding_exons_per_coding_transcript' => $$coding_stats[16],
                           'genebuild.stats.total_introns' => $$coding_stats[17],
                           'genebuild.stats.average_intron_length' => $$coding_stats[18]);

  my %noncoding_stats_hash = ('genebuild.stats.nc_non_coding_genes' => $$non_coding_stats[1],
                              'genebuild.stats.nc_small_non_coding_genes'  => $$non_coding_stats[2],
                              'genebuild.stats.nc_long_non_coding_genes' => $$non_coding_stats[3],
                              'genebuild.stats.nc_misc_non_coding_genes' => $$non_coding_stats[4],
                              'genebuild.stats.nc_average_genomic_span' => $$non_coding_stats[5],
                              'genebuild.stats.nc_average_sequence_length' => $$non_coding_stats[6],
                              'genebuild.stats.nc_shortest_gene_length' => $$non_coding_stats[7],
                              'genebuild.stats.nc_longest_gene_length' => $$non_coding_stats[8],
                              'genebuild.stats.nc_total_transcripts' => $$non_coding_stats[9],
                              'genebuild.stats.nc_transcripts_per_gene' => $$non_coding_stats[10],
                              'genebuild.stats.nc_total_exons' => $$non_coding_stats[11],
                              'genebuild.stats.nc_average_exon_length' => $$non_coding_stats[12],
                              'genebuild.stats.nc_average_exons_per_transcript' => $$non_coding_stats[13],
                              'genebuild.stats.nc_total_introns' => $$non_coding_stats[14],
                              'genebuild.stats.nc_average_intron_length' => $$non_coding_stats[15]);

  my %pseudogene_stats_hash = ('genebuild.stats.ps_pseudogenes' => $$pseudogene_stats[1],
                               'genebuild.stats.ps_average_genomic_span' => $$pseudogene_stats[2],
                               'genebuild.stats.ps_average_sequence_length' => $$pseudogene_stats[3],
                               'genebuild.stats.ps_shortest_gene_length' => $$pseudogene_stats[4],
                               'genebuild.stats.ps_longest_gene_length' => $$pseudogene_stats[5],
                               'genebuild.stats.ps_total_transcripts' => $$pseudogene_stats[6],
                               'genebuild.stats.ps_transcripts_per_gene' => $$pseudogene_stats[7],
                               'genebuild.stats.ps_total_exons' => $$pseudogene_stats[8],
                               'genebuild.stats.ps_average_exon_length' => $$pseudogene_stats[9],
                               'genebuild.stats.ps_average_exons_per_transcript' => $$pseudogene_stats[10],
                               'genebuild.stats.ps_total_introns' => $$pseudogene_stats[11],
                               'genebuild.stats.ps_average_intron_length' => $$pseudogene_stats[12]);

  my %assembly_stats_hash = ('assembly.stats.contig_n50' => $$assembly_stats[7],
                             'assembly.stats.total_genome_length' => $$assembly_stats[8],
                             'assembly.stats.total_coding_sequence_length' => $$assembly_stats[9],
                             'assembly.stats.total_gap_length' => $$assembly_stats[10],
                             'assembly.stats.spanned_gaps' => $$assembly_stats[11],
                             'assembly.stats.chromosomes' => $$assembly_stats[12],
                             'assembly.stats.toplevel_sequences' => $$assembly_stats[13],
                             'assembly.stats.component_sequences' => $$assembly_stats[14],
                             'assembly.stats.gc_percentage' => $$assembly_stats[15]);
  my $filename = $output_dir . "/stats_" . $dbname . ".sql";
  open(SQLOUT, ">", $filename)  or die "Cannot open file: $!";;
  print SQLOUT "USE " . $dbname .";\n"; 
  foreach(keys %coding_stats_hash) {
    if ($coding_stats_hash{$_} ne ''){
      print SQLOUT "INSERT INTO meta (species_id, meta_key, meta_value) VALUES(". $species_id .", '". $_ ."', ". $coding_stats_hash{$_} .");\n";
    }
  }
  foreach(keys %noncoding_stats_hash) {
    if ($noncoding_stats_hash{$_} ne ''){
      print SQLOUT "INSERT INTO meta (species_id, meta_key, meta_value) VALUES(". $species_id .", '". $_ ."', ". $noncoding_stats_hash{$_} .");\n";
    }
  }
  foreach(keys %pseudogene_stats_hash) {
    if ($pseudogene_stats_hash{$_} ne ''){
      print SQLOUT "INSERT INTO meta (species_id, meta_key, meta_value) VALUES(". $species_id .", '". $_ ."', ". $pseudogene_stats_hash{$_} .");\n";
    }
  }
foreach(keys %assembly_stats_hash) {
    if ($assembly_stats_hash{$_} ne ''){
      print SQLOUT "INSERT INTO meta (species_id, meta_key, meta_value) VALUES(". $species_id .", '". $_ ."', ". $assembly_stats_hash{$_} .");\n";
    }
  }
  close SQLOUT;
#}

my $coding_result = print_result(@final_coding_stats);
my $non_coding_result = print_result(@final_non_coding_stats);
my $pseudogene_result = print_result(@final_pseudogene_stats);
my $assembly_result = print_result(@final_assembly_stats);;

exit;


sub print_result {
  my (@results) = @_;

  my $results_string = "";
  foreach my $result (@results) {
    print Dumper($result);
    $results_string .= join("\t",@$result)."\n";
  }
  return($results_string);
}

sub process_coding_genes {
  my ($genes,$scientific_name) = @_;

  unless(scalar(@$genes)) {
    return([$scientific_name,0,'','','','','','','','','','','','','','','','','']);
  }

  my @sorted_results = ();

  my $species = $coding_headers[0];
  my $coding_genes_count = $coding_headers[1];
  my $avg_gene_genomic_length = $coding_headers[2];
  my $avg_gene_sequence_length = $coding_headers[3];
  my $avg_cds_length = $coding_headers[4];
  my $shortest_gene = $coding_headers[5];
  my $longest_gene = $coding_headers[6];
  my $transcript_count = $coding_headers[7];
  my $coding_transcript_count = $coding_headers[8];
  my $avg_transcripts_per_gene = $coding_headers[9];
  my $avg_coding_transcripts_per_gene = $coding_headers[10];
  my $exon_count = $coding_headers[11];
  my $coding_exon_count = $coding_headers[12];
  my $avg_exon_length = $coding_headers[13];
  my $avg_coding_exon_length = $coding_headers[14];
  my $avg_exons_per_transcript = $coding_headers[15];
  my $avg_coding_exons_per_coding_transcript = $coding_headers[16];
  my $intron_count = $coding_headers[17];
  my $avg_intron_length = $coding_headers[18];

  my $coding_sequence_length = 0;
  my $gene_count = scalar(@$genes);
  my $result = {$species => $scientific_name,
                $coding_genes_count => $gene_count,
                $avg_gene_genomic_length => 0,
                $avg_gene_sequence_length => 0,
                $avg_cds_length => 0,
                $shortest_gene => 999999999,
                $longest_gene => 0,
                $transcript_count => 0,
                $coding_transcript_count => 0,
                $avg_transcripts_per_gene => 0,
                $avg_coding_transcripts_per_gene => 0,
                $exon_count => 0,
                $coding_exon_count => 0,
                $avg_exon_length => 0,
                $avg_coding_exon_length => 0,
                $avg_exons_per_transcript => 0,
                $avg_coding_exons_per_coding_transcript => 0,
                $intron_count => 0,
                $avg_intron_length => 0,
               };

  my $total_genomic_length = 0;
  my $total_sequence_length = 0;
  my $total_cds_length = 0;
  my $total_coding_transcript_count = 0;
  my $total_transcript_count = 0;
  my $total_exons = 0;
  my $total_coding_exons = 0;
  my $total_introns = 0;
  my $total_exon_length = 0;
  my $total_coding_exon_length = 0;
  my $total_intron_length = 0;
  my $total_canonical_introns = 0;

  foreach my $gene (@$genes) {
    my $gene_length = $gene->length();
    if($gene_length < $result->{$shortest_gene}) {
      $result->{$shortest_gene} = $gene_length;
    }

    if($gene_length > $result->{$longest_gene}) {
      $result->{$longest_gene} = $gene_length;
    }

    $total_genomic_length += $gene_length;
    my $transcripts = $gene->get_all_Transcripts();
    my $canonical_transcript = $gene->canonical_transcript();
    $total_sequence_length += $canonical_transcript->length();
    foreach my $transcript (@$transcripts) {
      my $exons = $transcript->get_all_Exons();
      my $introns = $transcript->get_all_Introns();
      $total_exons += scalar(@$exons);
      $total_introns += scalar(@$introns);
      foreach my $exon (@$exons) {
        $total_exon_length += $exon->length();
      }

      foreach my $intron (@$introns) {
        $total_intron_length += $intron->length();
      }

      $total_transcript_count++;
      if($transcript->translation) {
       $total_coding_transcript_count++;
       # Only use the canonical for the adding to the coding seq length
        if($transcript->is_canonical()) {
          $coding_sequence_length += length($transcript->translateable_seq());
        }
        $total_cds_length += length($transcript->translateable_seq());
        my $coding_exons = $transcript->get_all_CDS();
        $total_coding_exons += scalar(@$coding_exons);
        foreach my $coding_exon (@$coding_exons) {
          $total_coding_exon_length += $coding_exon->length();
        } # End foreach my $coding_exon (@$coding_exons)
      } # if($transcript->translation)
    } # foreach my $transcript (@$transcripts)
  }

  $result->{$avg_gene_genomic_length} = sprintf("%.2f", $total_genomic_length/$gene_count);
  $result->{$avg_gene_sequence_length} = sprintf("%.2f", $total_sequence_length/$gene_count);
  $result->{$avg_cds_length} = sprintf("%.2f", $total_cds_length/$total_coding_transcript_count);
  $result->{$transcript_count} = $total_transcript_count;
  $result->{$coding_transcript_count} = $total_coding_transcript_count;
  $result->{$avg_transcripts_per_gene} = sprintf("%.2f", $total_transcript_count/$gene_count);
  $result->{$avg_coding_transcripts_per_gene} = sprintf("%.2f", $total_coding_transcript_count/$gene_count);
  $result->{$avg_exons_per_transcript} = sprintf("%.2f", $total_exons/$total_transcript_count);
  $result->{$exon_count} = $total_exons;
  $result->{$coding_exon_count} = $total_coding_exons;
  $result->{$avg_exon_length} = sprintf("%.2f", $total_exon_length/$total_exons);
  $result->{$avg_coding_exon_length} = sprintf("%.2f", $total_coding_exon_length/$total_coding_exons);
  $result->{$avg_exons_per_transcript} = sprintf("%.2f", $total_exons/$total_transcript_count);
  $result->{$avg_coding_exons_per_coding_transcript} = sprintf("%.2f", $total_coding_exons/$total_coding_transcript_count);
  $result->{$intron_count} = $total_introns;
  if($total_introns) {
    $result->{$avg_intron_length} = sprintf("%.2f", $total_intron_length/$total_introns);
  } else {
    $result->{$avg_intron_length} = '';
  }

  foreach my $key (keys(%$result)) {
    my $value = $result->{$key};
    my $sorted_index = first_index { $_ eq $key } @coding_headers;
    $sorted_results[$sorted_index] = $value;
  }

  return(\@sorted_results,$coding_sequence_length);
}


sub process_non_coding_genes {
  my ($genes,$scientific_name) = @_;

  unless(scalar(@$genes)) {
    return([$scientific_name,0,'','','','','','','','','','','','','','']);
  }

  my @sorted_results = ();

  my $species = $non_coding_headers[0];
  my $noncoding_genes_count = $non_coding_headers[1];
  my $snon_coding_count = $non_coding_headers[2];
  my $lnon_coding_count = $non_coding_headers[3];
  my $mnon_coding_count = $non_coding_headers[4];
  my $avg_gene_genomic_length = $non_coding_headers[5];
  my $avg_gene_sequence_length = $non_coding_headers[6];
  my $shortest_gene = $non_coding_headers[7];
  my $longest_gene = $non_coding_headers[8];
  my $transcript_count = $non_coding_headers[9];
  my $avg_transcripts_per_gene = $non_coding_headers[10];
  my $exon_count = $non_coding_headers[11];
  my $avg_exon_length = $non_coding_headers[12];
  my $avg_exons_per_transcript = $non_coding_headers[13];
  my $intron_count = $non_coding_headers[14];
  my $avg_intron_length = $non_coding_headers[15];

  my $gene_count = scalar(@$genes);
  my $result = {$species => $scientific_name,
                $noncoding_genes_count => $gene_count,
                $snon_coding_count => 0,
                $lnon_coding_count => 0,
                $mnon_coding_count => 0,
                $avg_gene_genomic_length => 0,
                $avg_gene_sequence_length => 0,
                $shortest_gene => 999999999,
                $longest_gene => 0,
                $transcript_count => 0,
                $avg_transcripts_per_gene => 0,
                $exon_count => 0,
                $avg_exon_length => 0,
                $avg_exons_per_transcript => 0,
                $intron_count => 0,
                $avg_intron_length => 0,
               };

  my $total_genomic_length = 0;
  my $total_sequence_length = 0;
  my $total_transcript_count = 0;
  my $total_exons = 0;
  my $total_introns = 0;
  my $total_exon_length = 0;
  my $total_intron_length = 0;
  my $total_canonical_introns = 0;

  foreach my $gene (@$genes) {
    my $gene_length = $gene->length();

    my $biotype_group = $gene->get_Biotype->biotype_group;
    if($biotype_group eq 'snoncoding') {
      $result->{$snon_coding_count}++;
    } elsif($biotype_group eq 'lnoncoding') {
      $result->{$lnon_coding_count}++;
    } elsif($biotype_group eq 'mnoncoding') {
      $result->{$mnon_coding_count}++;
    } else {
      die "Unexpected biotype group: ".$biotype_group;
    }

    if($gene_length < $result->{$shortest_gene}) {
      $result->{$shortest_gene} = $gene_length;
    }

    if($gene_length > $result->{$longest_gene}) {
      $result->{$longest_gene} = $gene_length;
    }

    $total_genomic_length += $gene_length;
    my $transcripts = $gene->get_all_Transcripts();
    my $canonical_transcript = $gene->canonical_transcript();
    $total_sequence_length += $canonical_transcript->length();
    foreach my $transcript (@$transcripts) {
      my $exons = $transcript->get_all_Exons();
      my $introns = $transcript->get_all_Introns();
      $total_exons += scalar(@$exons);
      $total_introns += scalar(@$introns);
      foreach my $exon (@$exons) {
        $total_exon_length += $exon->length();
      }

      foreach my $intron (@$introns) {
        $total_intron_length += $intron->length();
      }

      $total_transcript_count++;
    } # foreach my $transcript (@$transcripts)
  }

  $result->{$avg_gene_genomic_length} = sprintf("%.2f", $total_genomic_length/$gene_count);
  $result->{$avg_gene_sequence_length} = sprintf("%.2f", $total_sequence_length/$gene_count);
  $result->{$transcript_count} = $total_transcript_count;
  $result->{$avg_transcripts_per_gene} = sprintf("%.2f", $total_transcript_count/$gene_count);
  $result->{$avg_exons_per_transcript} = sprintf("%.2f", $total_exons/$total_transcript_count);
  $result->{$exon_count} = $total_exons;
  $result->{$avg_exon_length} = sprintf("%.2f", $total_exon_length/$total_exons);
  $result->{$avg_exons_per_transcript} = sprintf("%.2f", $total_exons/$total_transcript_count);
  $result->{$intron_count} = $total_introns;
  if($total_introns) {
    $result->{$avg_intron_length} = sprintf("%.2f", $total_intron_length/$total_introns);
  } else {
    $result->{$avg_intron_length} = '';
  }

  foreach my $key (keys(%$result)) {
    my $value = $result->{$key};
    my $sorted_index = first_index { $_ eq $key } @non_coding_headers;
    $sorted_results[$sorted_index] = $value;
  }

  return(\@sorted_results);
}

sub process_pseudogenes {
  my ($genes,$scientific_name) = @_;

  unless(scalar(@$genes)) {
    return([$scientific_name,0,'','','','','','','','','','','']);
  }

  my @sorted_results = ();

  my $species = $pseudogene_headers[0];
  my $pseudogenes_count = $pseudogene_headers[1];
  my $avg_gene_genomic_length = $pseudogene_headers[2];
  my $avg_gene_sequence_length = $pseudogene_headers[3];
  my $shortest_gene = $pseudogene_headers[4];
  my $longest_gene = $pseudogene_headers[5];
  my $transcript_count = $pseudogene_headers[6];
  my $avg_transcripts_per_gene = $pseudogene_headers[7];
  my $exon_count = $pseudogene_headers[8];
  my $avg_exon_length = $pseudogene_headers[9];
  my $avg_exons_per_transcript = $pseudogene_headers[10];
  my $intron_count = $pseudogene_headers[11];
  my $avg_intron_length = $pseudogene_headers[12];

  my $gene_count = scalar(@$genes);
  my $result = {$species => $scientific_name,
                $pseudogenes_count => $gene_count,
                $avg_gene_genomic_length => 0,
                $avg_gene_sequence_length => 0,
                $shortest_gene => 999999999,
                $longest_gene => 0,
                $transcript_count => 0,
                $avg_transcripts_per_gene => 0,
                $exon_count => 0,
                $avg_exon_length => 0,
                $avg_exons_per_transcript => 0,
                $intron_count => 0,
                $avg_intron_length => 0,
               };

  my $total_genomic_length = 0;
  my $total_sequence_length = 0;
  my $total_transcript_count = 0;
  my $total_exons = 0;
  my $total_introns = 0;
  my $total_exon_length = 0;
  my $total_intron_length = 0;
  my $total_canonical_introns = 0;

  foreach my $gene (@$genes) {
    my $gene_length = $gene->length();

    if($gene_length < $result->{$shortest_gene}) {
      $result->{$shortest_gene} = $gene_length;
    }

    if($gene_length > $result->{$longest_gene}) {
      $result->{$longest_gene} = $gene_length;
    }

    $total_genomic_length += $gene_length;
    my $transcripts = $gene->get_all_Transcripts();
    my $canonical_transcript = $gene->canonical_transcript();
    $total_sequence_length += $canonical_transcript->length();
    foreach my $transcript (@$transcripts) {
      my $exons = $transcript->get_all_Exons();
      my $introns = $transcript->get_all_Introns();
      $total_exons += scalar(@$exons);
      $total_introns += scalar(@$introns);
      foreach my $exon (@$exons) {
        $total_exon_length += $exon->length();
      }

      foreach my $intron (@$introns) {
        $total_intron_length += $intron->length();
      }

      $total_transcript_count++;
    } # foreach my $transcript (@$transcripts)
  }

  $result->{$avg_gene_genomic_length} = sprintf("%.2f", $total_genomic_length/$gene_count);
  $result->{$avg_gene_sequence_length} = sprintf("%.2f", $total_sequence_length/$gene_count);
  $result->{$transcript_count} = $total_transcript_count;
  $result->{$avg_transcripts_per_gene} = sprintf("%.2f", $total_transcript_count/$gene_count);
  $result->{$avg_exons_per_transcript} = sprintf("%.2f", $total_exons/$total_transcript_count);
  $result->{$exon_count} = $total_exons;
  $result->{$avg_exon_length} = sprintf("%.2f", $total_exon_length/$total_exons);
  $result->{$avg_exons_per_transcript} = sprintf("%.2f", $total_exons/$total_transcript_count);
  $result->{$intron_count} = $total_introns;
  if($total_introns) {
    $result->{$avg_intron_length} = sprintf("%.2f", $total_intron_length/$total_introns);
  } else {
    $result->{$avg_intron_length} = '';
  }

  foreach my $key (keys(%$result)) {
    my $value = $result->{$key};
    my $sorted_index = first_index { $_ eq $key } @pseudogene_headers;
    $sorted_results[$sorted_index] = $value;
  }

  return(\@sorted_results);
}
