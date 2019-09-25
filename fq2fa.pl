use strict;
use warnings;
use File::Basename;
use Cwd qw(abs_path);

die "perl $0 <in.fq> <out.fa>\n" if (@ARGV != 2);

fq2fa();

sub fq2fa
{
  open TM, "> $ARGV[1].tmp" or die $!;
  my ($read, $flag) = ('', 0);
  open IN, ($ARGV[0] =~ /\.gz$/) ? "gunzip -dc $ARGV[0] | " : $ARGV[0] or die $!;
  while(<IN>)
  {
    chomp;
    if ($flag == 0) {
      $read = $_;
      $read = $1 if (/.*\_U(\d+).*/);
    }
    print TM "$_\t$read\n" if ($flag == 1);
    $flag ++;
    $flag = 0 if ($flag == 4);
  }
  close IN;
  close TM;
  my ($pre, $num, $dir) = ('na', 0, dirname(abs_path($ARGV[1])));
  $flag = 1;
  my %id; #--- count UMI
  open IT, "sort --buffer-size=1000M --temporary-directory=$dir -k 1,1 $ARGV[1].tmp | " or die $!;
  open OT, "> $ARGV[1]" or die $!;
  while (<IT>) 
  {
    chomp;
    my @tmp = split /\t/, $_;
    if ($tmp[0] ne $pre) {
      if ($pre ne 'na') {
        print OT ">tag$flag\_$num\_", scalar keys %id, "\n$pre\n";
        $flag ++;
      }
      ($pre, $num) = ($tmp[0], 0);
      %id = ();
    }
    $num ++;
    $id{$tmp[1]} = '';
  }
  print OT ">tag$flag\_$num\_", scalar keys %id, "\n$pre";
  close OT;
  close IT;
  unlink "$ARGV[1].tmp";
}

exit 1;
