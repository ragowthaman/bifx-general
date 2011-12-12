#!/depot/perl-5.12.1/bin/perl
#
#
# PURPOSE: Take a multi fasta, blast each sequence aginst a database and 
#          render a graphical display for each blast search 
#          Uses NCBI's blastp/blastn output
#
#		   Also takes two single fasta and compare& render using bl2seq
#		   (set -bl2seq & supply -j)
# 
#		   Also renders a done single query blastoutfile. (supply -blastfile)
#
# Author: gowthaman.ramasamy@seattlebiomed.org
# Core of the code is liberally stollen from Lincoln Stein's :http://www.bioperl.org/wiki/HOWTO:Graphics
#
################################################################################################################
use strict;
use warnings;
use DateTime;
use Bio::SeqIO;
use Math::Round;
use Data::Dumper;
use Bio::Graphics;
use Bio::SearchIO;
use Bio::Tools::GFF;
#use Bio::Tools::Run::StandAloneBlast;
use Bio::SeqFeature::Generic;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);


#----------------------------------------------------------------------
# handle command line
#----------------------------------------------------------------------
my $usage =qq{
$0 [options] -i <fasta file name> -d <database> -p <blast program>


OPTIONS:
    -i, --input
        Fasta file to be quried (same as blastall's -i)
    -d, --database
        Blast database name (same as blastall's -d)
    -p, --program
    	Blastprogram (same as blastall's -p)
    --prefix
        prefix for the png outfile name
    -o, --outformat
        format of the grapical outfile. Default: png
    -e, --evalue
        Evalue cutoff (same as blastall's -e)
    
    BL2SEQ options:
    --bl2seq
    	If bl2seq needs to be run
    -j, --input2 in case of bl2seq
        Fasta file to be used as subject in bl2seq. Ingnore for other blast.
    
    RENDERING BLASTOUTFILE
    -blastfile
    	Will render blastfile directly. NO blast is done.

    -help
        this help message
    -debug
};
my %options =();
my $opt_obj = GetOptions (\%options,
              'input|i=s',
              'database|d=s',
              'program|p=s',
              'prefix=s',
              'outformat|o=s',
              'evalue|e=s',
              'bl2seq',
              'input2|j=s',
              'blastfile=s',
   			  'debug|d',
              'help|h',
              );

die("$usage\n") if($options{help});

#----------------------------------------------------------------------
# Handle command line
#----------------------------------------------------------------------
$options{outformat} ||= "png";
$options{prefix} ||= 'analysis';

my $dt = DateTime->now();
my $today = $dt->ymd;
#----------------------------------------------------------------------
# Main - main
#----------------------------------------------------------------------

print Dumper(\%options) if $options{debug};
#read and iterate thru a multifasta
# walk thru target fasta
my $in = Bio::SeqIO->new(-format    => 'fasta',
                         -file      => $options{input}) || die("Could not open file $options{fasta}");

if($options{blastfile}){
	my $graphicfile = $options{blastfile}."_graph.".$options{outformat};
	&render_blast($options{blastfile}, $graphicfile, $options{outformat});
	
}elsif($options{bl2seq}){
	my $blastfile = $options{input}.'_'.$options{input2}."_blast.out";
	my $graphicfile = $blastfile."_graph.".$options{outformat};
    
	&bl2seq($options{input}, $options{input2}, $blastfile);
	&render_blast($blastfile, $graphicfile, $options{outformat}); 
	
}else{
	while (my $seq = $in->next_seq){
		my $qid = $seq->id;
		my $sequence = $seq->seq();
		
		my $seqfile = $qid."_seq.fasta";
		my $blastfile = $options{prefix}.'_'.$qid."_blast.out";
		my $graphicfile = $options{prefix}.'_'.$qid."_graph.".$options{outformat};
		
		# write individual fasta
		system("rm $seqfile") if -e $seqfile;
		open(my $FH, ">>$seqfile");
		print $FH ">$qid\n$sequence\n";
		close($FH);
		
		#blast
		&blast($seqfile, $blastfile);
		
		#render graphics
		&render_blast($blastfile, $graphicfile, $options{outformat});    
	}
}

sub blast(){
	my ($infile, $blastfile) = @_;
	my $cmd =qq{blastall -p  $options{program} -d $options{database} -i $infile -o $blastfile};
	print $cmd if $options{debug};
	$cmd .= " -e $options{evalue} " if $options{evalue};
	system("$cmd");
}

sub bl2seq(){
	my ($infile, $infile2, $blastfile) = @_;
	my $cmd =qq{bl2seq -p  $options{program} -i $infile -j $infile2 -o $blastfile};
	$cmd .= " -e $options{evalue} " if $options{evalue};
	system("$cmd");
}


sub render_blast(){
	my ($blastfile, $outfile, $outformat) = @_;
	
	my $searchio = Bio::SearchIO->new(-file   => $blastfile,
                                  -format => 'blast') or die "parse failed";
 
my $result = $searchio->next_result() or die "no result";
 
my $panel = Bio::Graphics::Panel->new(
                                      -length    => $result->query_length,
                                      -width     => 800,
                                      -pad_left  => 10,
                                      -pad_right => 10,
                                     );
 
my $full_length = Bio::SeqFeature::Generic->new(
                                                -start        => 1,
                                                -end          => $result->query_length,
                                                -display_name => $result->query_name,
                                               );
$panel->add_track($full_length,
                  -glyph   => 'arrow',
                  -tick    => 2,
                  -fgcolor => 'black',
                  -double  => 1,
                  -label   => 1,
                 );
 
my $track = $panel->add_track(
                              -glyph       => 'graded_segments',
                              -label       => 1,
                              -connector   => 'dashed',
                              -bgcolor     => 'blue',
                              -font2color  => 'red',
                              -sort_order  => 'high_score',
                              -description => sub {
                                my $feature = shift;
                                return unless $feature->has_tag('description');
                                my ($description) = $feature->each_tag_value('description');
                                my $score = $feature->score;
                                my ($length) = $feature->each_tag_value('hit_len');
                                my ($tot_hsplen) = $feature->each_tag_value('tot_hsplen');
                                my $tot_hspcov =round(($tot_hsplen/$length)*100);
                                $tot_hspcov = $tot_hspcov.'%';
                                

                                "$description, score=$score, Sub_tot_length=$length, tot_hsplen=$tot_hsplen, tot_hspcov=$tot_hspcov";
                               },
                             );
 
while( my $hit = $result->next_hit ) {
  #next unless $hit->significance < 1E-20;
  my $feature = Bio::SeqFeature::Generic->new(
                                              -score        => $hit->raw_score,
                                              -display_name => $hit->name,
                                              -tag          => {
                                                                description => $hit->description,
                                                                hit_len => $hit->length,
                                                               },
                                              
                                             );
  my $total_hsp_length = "0";
  while( my $hsp = $hit->next_hsp ) {
  	 
    $feature->add_sub_SeqFeature($hsp,'EXPAND');
    $total_hsp_length += $hsp->length('hit');
  }
  $feature->add_tag_value('tot_hsplen',$total_hsp_length);

  $track->add_feature($feature);
}
open(FH, ">>$outfile"); 
print FH $panel->png;
close(FH);
}