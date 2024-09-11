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

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> July 2021

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
		    outversion => 'N0016',
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

mkpath($opts{outdir});

sub mod_N0012_hv1_update {
  my ($wav, $qe) = @_;

=begin comment

From Nick, 2012-10-30:

Hi Pete,
This is the correction line for the HRC QEU [sic] as a function of wavelength.

1/correction = 1.04344 + (9.40064e-5 * lambda)

You should end up with a 5.9% decrease in QE at 165AA and 4.4% @ 10AA.
This is effective for all hrcs after the voltage change, from March 29 -
on.

=cut

  my ($a, $b) = (9.40064e-5, 1.04344);
  $qe *= $a * $wav + $b;

}

sub mod_N0013_hi_e_corr {
  my ($wav, $qe) = @_;

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

  for my $i (0..$wav->getdim(1)-1) {

    # FIXME: we want a 1.07 raise in flux, so for the next iteration
    # the factor should be 1/1.07 instead of 0.93
    # my ($factor, $w1, $w2) = (0.93, 33.9, 45);
    my ($factor, $w1, $w2) = (1/1.07, 33.9, 45);

    my $wav = $wav->slice(",($i)");

    my $corr = ones($wav->nelem);

    my $index = which($wav<$w1);
    (my $tmp = $corr->index($index)) .= $factor;

    $index = which(($wav>=$w1) & ($wav<$w2));
    ($tmp = $corr->index($index)) .= $factor + (1-$factor)/($w2-$w1) * ($wav->index($index) - $w1);

    ($tmp = $qe->slice(",($i)")) *= $corr;
  }
}

sub mod_N0014_lsfparm_update {
  my ($wav, $qe) = @_;

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
  $qe /= interpol($wav, $ratio_wav, $ratio);
}

sub mod_N0016_hv2_update {
  my $qe = shift;
  $qe *= 1.106;
}

my @cvsd = qw/ 1999-07-22 2012-03-29 2021-05-14 /;

for my $j (0..$#cvsd) {

  my $infile = $opts{infile};

  my $cvsd = $cvsd[$j];

  my $outfile = "$opts{outdir}/hrcsD${cvsd}qe$opts{outversion}.fits";

  print STDERR "cp $infile -> $outfile\n";
  copy($infile, $outfile) or die $!;

  my $s = 0;
  my $fptr = Astro::FITS::CFITSIO::open_file($outfile, Astro::FITS::CFITSIO::READWRITE(), $s);

  update_keywords($fptr, $cvsd, $s);

  $fptr->write_date($s);
  $fptr->write_chksum($s);

  for my $hdrnum (2..4) {

    $fptr->movabs_hdu($hdrnum, Astro::FITS::CFITSIO::BINARY_TBL(), $s);
    update_keywords($fptr, $cvsd, $s);

    my ($e, $qe) = read_bintbl_cols($fptr, qw( energy qe ));
    my $wav = 12.39854 / $e;

    mod_N0012_hv1_update($wav, $qe) if $cvsd gt '1999-07-22';
    mod_N0013_hi_e_corr($wav, $qe);
    mod_N0014_lsfparm_update($wav, $qe);
    mod_N0016_hv2_update($qe) if $cvsd ge '2021-05-14';

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
  my ($fptr, $cvsd, $s) = @_;

  $fptr->update_key_str('creator', $FindBin::Bin . '/' . $FindBin::RealScript, undef, $s);

  $fptr->update_key_str('cvsd0001', $cvsd.'T00:00:00', undef, $s);
  $fptr->update_key_str('cvst0001', '00:00:00', undef, $s);

  check_status($s);

}
