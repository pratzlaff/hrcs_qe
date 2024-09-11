#! /usr/bin/perl

use strict;
use warnings;

=head1 NAME

template - A template for Perl programs.

=head1 SYNOPSIS

cp template newprog

=head1 DESCRIPTION

blah blah blah

=head1 OPTIONS

=over 4

=item --help

Show help and exit.

=item --version

Show version and exit.

=back

=head1 AUTHOR

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> Mar 2020

=head1 SEE ALSO

perl(1).

=cut

my $version = '0.1';

use FindBin;
use Config;
use Carp;
use Astro::FITS::CFITSIO qw( CASEINSEN );
use Chandra::Tools::Common qw( read_bintbl_cols );
use File::Copy;
use File::Path qw/ mkpath /;
use PDL;
use PDL::Graphics::PGPLOT;
use PDL::IO::Misc;
use PDL::NiceSlice;

use Getopt::Long;
my %default_opts = (
		    infile => '../N0011/hrcsD1999-07-22qeN0011.fits',
		    outversion => 'N0015',
		    outdir => './qe',
		    ratiofile => '../N0014/NewOldRatio.out',
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'infile=s', 'outfile=s',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

my $infile = $opts{infile};

mkpath($opts{outdir});

# for years later than 2012.5, the maximum correction is 2.1%/yr at 0
# Ang, and everything from 20 AA down to 0 AA linearly transisitions
# from 1.0 to that maximum
sub short_corr {
  my ($y, $wav) = @_;
  my $y0 = 2012.5;
  my ($w0, $w1) = (0, 20);
  my $rate = 0.021;

  my $corr = ones($wav);
  return $corr if $y < $y0;

  my $max = 1 - ($y-$y0) * $rate;

  my $tmp;
  my $ii = which(($wav>0) & ($wav<$w0));
  ($tmp = $corr->index($ii)) .= $max;
  $ii = which(($wav>=$w0) & ($wav<$w1));
  ($tmp = $corr->index($ii)) .= 1 - (1-$max) * ($w1-$wav->index($ii)) / ($w1-$w0);

  return $corr;
}

=begin comment

From Brad's 2015-06-26 email:
-----------------------------

The ratio newEEFRAC/oldEEFRAC for the standard bowtie region is in
/data/letg4/bradw/EEFRACS2013/EEFRAC/LSFPARMstuff/NewOldRatio.out
col1 = |wavelength|
col2 = new/old for - wavelengths
col3 = new/old for + wavelengths
(Ignore cols 4 and 5--was for debugging)

=cut

my ($ratio_wav, $ratio_neg, $ratio_pos) = rcols($opts{ratiofile}, 0..2, { lines => '2:-1' });

# the $ratio_neg and $ratio_pos are typically within 0.5%, so just use
# the average to make everything easier
my $ratio = 0.5 * ($ratio_neg + $ratio_pos);

# technically should be something like 2012.625 (half way between 2012-03-29 and 2013-01-01),
# but this small correction isn't justified by the mkn421 w3-6_E4 plot
my @cvsd = qw/ 1999-07-22 2012-03-29 /;
my @years = ( 1999.56, 2012.24 );

# these will have the short wavelength correction on S2
for my $y (2013..2022) {
  push @cvsd, "$y-01-01";
  push @years, $y+0.5;
}

for my $j (0..$#cvsd) {

  my $cvsd = $cvsd[$j];
  my $year = $years[$j];

  my $cved;
  $cved = $cvsd[$j+1] unless $j == $#cvsd;

  my $outfile = "$opts{outdir}/hrcsD${cvsd}qe$opts{outversion}.fits";

  print STDERR "cp $infile -> $outfile\n";
  copy($infile, $outfile) or die $!;

  my $s = 0;
  my $fptr = Astro::FITS::CFITSIO::open_file($outfile, Astro::FITS::CFITSIO::READWRITE(), $s);


# Key   72: C              CVSD0001     = 1999-07-22T00:00:00  /
# Key   73: C              CVST0001     = 00:00:00             /


  update_keywords($fptr, $cvsd, $cved, $s);

  $fptr->write_date($s);
  $fptr->write_chksum($s);

  for my $hdrnum (2..4) {

    $fptr->movabs_hdu($hdrnum, Astro::FITS::CFITSIO::BINARY_TBL(), $s);

    update_keywords($fptr, $cvsd, $cved, $s);

    my ($energy, $qe) = read_bintbl_cols($fptr, qw( energy qe ));

    my $wavelength = 12.39854 / $energy;

=begin comment

From Jeremy's 2014-11-25 email:

Hello Pete, for the high E end QE correction, the recipe for the last
correction is on

http://cxc.harvard.edu/twiki/bin/view/HrcCal/Letg_lo

look for "QE Correction Curve"

the summary from there is:

    No change below 33.9AA
    Linear interpolation from 33.9 - 45AA
    Based on PKS2155-304 45 - 70AA
    Based on HZ43 > 70AA
    Removed spurious feature near 168AA

So, we need to make the change below 33.9 AA where there was no change
before, and linearly interpolate to 1 (ie no correction) in the same
way from 33.9 to 45.

I I recall, (though check the sign!), we needed a 7% (or whatever the
number was) increase in the HRC-S flux to match the ACIS-S flux.  This
would mean a 7% reduction in QE, which (if the sign is correct!)
meshes well with our earlier correction.  Earlier, we reduced the QE
at longer wavelengths, but kept shorter the same. Now we see we need
to reduce at shorter wavelengths too.  This is good (bad if I have the
sign wrong) because it means we don;t (do) introduce larger QE
gradients compared with earlier times.

So, hopefully our QE correction factor will be 0.93 from 1-33.9AA and
ramping up to 1 at 45.

=cut

    for my $i (0..$energy->getdim(1)-1) {
      # FIXME: we want a 1.07 raise in flux, so for the next iteration
      # the factor should be 1/1.07 instead of 0.93
#      my ($factor, $w1, $w2) = (0.93, 33.9, 45);
      my ($factor, $w1, $w2) = (1/1.07, 33.9, 45);

      my $wavelength = $wavelength->slice(",($i)");

      my $corr = ones($wavelength->nelem);

      my $index = which($wavelength<$w1);
      (my $tmp = $corr->index($index)) .= $factor;

      $index = which(($wavelength>=$w1) & ($wavelength<$w2));
      ($tmp = $corr->index($index)) .= $factor + (1-$factor)/($w2-$w1) * ($wavelength->index($index) - $w1);

      ($tmp = $qe->slice(",($i)")) *= $corr;

      if ($hdrnum == 3) {
	($tmp = $qe->slice(",($i)")) /= short_corr($year, $wavelength);
      }

    }

=begin comment

From Nick, 2012-10-30:

Hi Pete,
This is the correction line for the HRC QEU [sic] as a function of wavelength.

1/correction = 1.04344 + (9.40064e-5 * lambda)

You should end up with a 5.9% decrease in QE at 165AA and 4.4% @ 10AA.
This is effective for all hrcs after the voltage change, from March 29 -
on.

=cut

    if ($cvsd ne '1999-07-22') {
      my ($a, $b) = (9.40064e-5, 1.04344);

      $qe *= ($a * $wavelength + $b);
    }

    $qe /= interpol($wavelength, $ratio_wav, $ratio);

    my $qe_colnum;
    $fptr->get_colnum(CASEINSEN, 'qe', $qe_colnum, $s);
    $fptr->write_col(Astro::FITS::CFITSIO::TDOUBLE(), $qe_colnum, 1, 1, $qe->nelem, $qe->double->get_dataref, $s);

    $fptr->write_date($s);
    $fptr->write_chksum($s);

    check_status($s);

  }

  $fptr->close_file($s);
}

exit 0;

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}

sub check_status {
  my $s = shift;
  if ($s != 0) {
    my $txt;
    Astro::FITS::CFITSIO::fits_get_errstatus($s,$txt);
    carp "CFITSIO error: $txt";
    return 0;
  }

  return 1;
}

sub update_keywords {
  my ($fptr, $cvsd, $cved, $s) = @_;

  $fptr->update_key_str('creator', $FindBin::Bin . '/' . $FindBin::RealScript, undef, $s);

  $fptr->update_key_str('cvsd0001', $cvsd.'T00:00:00', undef, $s);
  $fptr->update_key_str('cvst0001', '00:00:00', undef, $s);

  if (defined $cved) {
    $fptr->update_key_str('cved0001', $cved, undef, $s);
    $fptr->update_key_str('cvet0001', '00:00:00', undef, $s);
  }

  check_status($s);

}
