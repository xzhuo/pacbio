use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

GetOptions(
    "code"      => \ my $code,  # use CodecV1 instead of frame for calculation. Default option is to convert the code to frame number.
    "sam=s"     => \ my $sam,
    "window=i"  => \(my $win = 10),
    "out=s"     => \ my $out
) or die "Invalid options passed to $0\n";


if($sam=~/\.bam$/) {
    open IN, "samtools view -\@5 $sam|";
}
elsif($sam=~/\.sam$/) {
    open IN, $sam;
}
else {
    die "$sam is no a bam/sam file!\n";
}

while(<IN>)
{
    chomp;
    next if /^@/;
    my $line=$_;
    my @F=split /\t/, $_, 12;
    my ($query, $flag, $ref, $mapPos, $mapQ, $cigar, $rnext, $pnext, $tlen, $seq, $qual) = @F;
    my %tags=split /:\S:|\s+/, $F[11];

    my $rg = $tags{"RG"};
    my $ec = $tags{"ec"};  #  Effective coverage
    my $fn = $tags{"fn"};  # Forward number of complete passes
    my $np = $tags{"np"};  # number of passes.
    my $rn = $tags{"rn"};  # Reverse number of complete passes
    my $rq = $tags{"rq"};  # Float in [0, 1] encoding expected accuracy
    my $zm = $tags{"zm"};
    my @fi = split /,/, $tags{"fi"};
    my @fp = split /,/, $tags{"fp"};
    my @ri = split /,/, $tags{"ri"};
    my @rp = split /,/, $tags{"rp"};
    my @sn = split /,/, $tags{"sn"};  # 4 floats for the average signal-to-noise ratio of A, C, G, and T (in that order) over the HQRegion
    shift @fi;
    shift @fp;
    shift @ri;
    shift @rp;
    unless ($code) {
        @fi = map {load_frame_codes($_)} @fi;
        @fp = map {load_frame_codes($_)} @fp;
        @ri = map {load_frame_codes($_)} @ri;
        @rp = map {load_frame_codes($_)} @rp;
    }
    @ri = reverse @ri;
    @rp = reverse @rp;
    shift @sn; 
    my @cpgSites = indexes($seq, "CG");
    # print Dumper(@cpgSites);
    foreach my $i(@cpgSites) {
        my $start = $i - $win;
        my $end = $i + $win + 1;
        my $length = $end - $start + 1;
        next if $start < 0 || $end >= length($seq);
        my $substring = substr($seq, $start, $length);
        my $subfi = join ",", @fi[$start..$end];
        my $subfp = join ",", @fp[$start..$end];
        my $subri = join ",", @ri[$start..$end];
        my $subrp = join ",", @rp[$start..$end];
        # print Dumper($substring);
        print join("\t", $query, $start, $end, "CG", $length, $np, $substring, $subfi, $subfp, $subri, $subrp), "\n";
    }
}
close IN;


sub indexes {
    my $string = shift;
    my $sub = shift;
    my $offset = 0;
    my @out = ();
    my $index = index($string, $sub, $offset);
    while ($index != -1) {
        push @out, $index;
        $offset = $index + 1;
        $index = index($string, $sub, $offset);
    }
    return @out;
}

sub load_frame_codes {
    my $code = shift;
    my $frame;
    if ($code < 64) {
        $frame = $code;
    }
    elsif ($code < 128) {
        $frame = 64 + ($code - 64)*2;
    }
    elsif ($code < 192) {
        $frame = 192 + ($code - 128)*4;
    }
    elsif ($code < 255) {
        $frame = 448 + ($code - 192)*8;
    }
    else {
        die "ip or pw frame code wrong";
    }
    return $frame;
}

