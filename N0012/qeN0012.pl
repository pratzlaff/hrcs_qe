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

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> November 2012

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
use PDL;
use PDL::Graphics::PGPLOT;
use PDL::IO::Misc;
use PDL::NiceSlice;

use Getopt::Long;
my %default_opts = (
		    infile => 'hrcsD1999-07-22qeN0011.fits',
		    outfile1 => 'hrcsD1999-07-22qeN0012.fits',
		    outfile2 => 'hrcsD2012-03-29qeN0012.fits',
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

my ($infile, $outfile1, $outfile2) = @opts{qw( infile outfile1 outfile2 )};

for my $outfile ($outfile1, $outfile2) {

  print STDERR "cp $infile -> $outfile\n";
  copy($infile, $outfile) or die $!;

  my $s = 0;
  my $fptr = Astro::FITS::CFITSIO::open_file($outfile, Astro::FITS::CFITSIO::READWRITE(), $s);


# Key   72: C              CVSD0001     = 1999-07-22T00:00:00  /
# Key   73: C              CVST0001     = 00:00:00             /

  my ($cvsd, $cved, $cvst, $cvet);
  if ($outfile eq $outfile1) {
    $cvsd = sprintf '%04d-%02d-%02dT00:00:00', 1999, 7, 22;
    $cvst = '00:00:00';

    $cved = sprintf '%04d-%02d-%02dT23:59:59', 2012, 3, 29;
    $cvet = '00:00:00';
  }
  else {
    $cvsd = sprintf '%04d-%02d-%02dT00:00:00', 2012, 3, 29;
    $cvst = '00:00:00';
  }

  update_keywords($fptr, $cvsd, $cvst, $cved, $cvet, $s);

  $fptr->write_date($s);
  $fptr->write_chksum($s);

  for my $hdrnum (2..4) {

    $fptr->movabs_hdu($hdrnum, Astro::FITS::CFITSIO::BINARY_TBL(), $s);

    update_keywords($fptr, $cvsd, $cvst, $cved, $cvet, $s);

=begin comment

From Nick, 2012-10-30:

Hi Pete,
This is the correction line for the HRC QEU as a function of wavelength.

1/correction = 1.04344 + (9.40064e-5 * lambda)

You should end up with a 5.9% decrease in QE at 165AA and 4.4% @ 10AA.
This is effective for all hrcs after the voltage change, from March 29 -
on.

=cut

    if ($outfile eq $outfile2) {
      my $hdr = $fptr->read_header();
      my ($a, $b) = (9.40064e-5, 1.04344);

      my ($energy, $qe) = read_bintbl_cols($fptr, qw( energy qe ));
      my $lam = 12.39854 / $energy;

      $qe *= ($a * $lam + $b);

      my $qe_colnum;
      $fptr->get_colnum(CASEINSEN, 'qe', $qe_colnum, $s);
      $fptr->write_col(Astro::FITS::CFITSIO::TDOUBLE(), $qe_colnum, 1, 1, $qe->nelem, $qe->double->get_dataref, $s);

    }

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
  my ($fptr, $cvsd, $cvst, $cved, $cvet, $s) = @_;

  $fptr->update_key_str('creator', $FindBin::Bin . '/' . $FindBin::RealScript, undef, $s);

  $fptr->update_key_str('cvsd0001', $cvsd, undef, $s);
  $fptr->update_key_str('cvst0001', $cvst, undef, $s);

  if (defined $cved) {
    $fptr->update_key_str('cved0001', $cved, undef, $s);
    $fptr->update_key_str('cvet0001', $cvet, undef, $s);
  }

  check_status($s);

}
