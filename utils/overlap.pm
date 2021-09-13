        package overlap;

        sub overlap {
            my $startlocus1  = $_[0];
            my $endlocus1    = $_[1];
            my $startlocus2  = $_[2];
            my $endlocus2    = $_[3];
            my $locus1chr    = $_[4];
            my $locus2chr    = $_[5];
            my $locus1strand = $_[6];
            my $locus2strand = $_[7];
            my $overlap;
            die "overlap::overlap subroutine expects 8 arguments\n"
                unless ( $#_ == 7 );
            die "overlap::overlap subroutine: strand must be '+' or '-':"
                . join( " ", @_ )
                . "Died.\n"
                unless ( ( $locus1strand eq '+' || $locus1strand eq '-' )
                && ( $locus2strand eq '+' || $locus2strand eq '-' ) );

            $locus1strand = '+' if ( $locus1strand eq '.' );
            $locus2strand = '+' if ( $locus2strand eq '.' );

# THIS SCRIPT PERFORMS STRANDED COMPARISONS ONLY. if 2 transcripts overlap but are on different strands, return value will be 0.
# "up/downstream" is to be taken in the transcriptional sense.
# It returns a 3-element array where:
# $returnedArray[0] is the overlap class (see below)
#	$returnedArray[1] is the distance between 5' ends of the two loci being compared.
#           < 0 if locus1's 5' end is upstream of locus2's
#           > 0 if locus1's 5' end is downstream of locus2's
#           undef if $returnedArray[0] == 0
# $returnedArray[2] is the distance between 3' ends of the two loci being compared.
#           < 0 if locus1's 3' end is upstream of locus2's
#           > 0 if locus1's 3' end is downstream of locus2's
#           undef if $returnedArray[0] == 0

            # $returnedArray[0] return codes:
            #
            ###    0           if locus1 and 2 are on different chromosomes (no overlap) and/or different strands

            ###    -1          if locus 2 is upstream (5') of locus 1 (no overlap)

          #               locus2                 >>>>>>>>>
          #                      -------------------------------------
          #               locus1                              >>>>>>>>
          # or ===============================================================
          #               locus2                 <<<<<<<<<
          #                      -------------------------------------
          #               locus1       <<<<<<<<

            ###    -2          if locus 2 is downstream (3') of locus 1 (no overlap)

          #               locus2                 >>>>>>>>>
          #                      -------------------------------------
          #               locus1       >>>>>>>>
          # or ===============================================================
          #               locus2                 <<<<<<<<<
          #                      -------------------------------------
          #               locus1                              <<<<<<<<

            ###    1           if locus1 is fully included in locus2

          #               locus2                 >>>>>>>>>
          #                      -------------------------------------
          #               locus1                   >>>>
          # or ===============================================================
          #               locus2                 <<<<<<<<<
          #                      -------------------------------------
          #               locus1                   <<<<

            ###    2           if locus1 3' is inside locus2, but locus1 5' extends locus2 5'

          #               locus2                 >>>>>>>>>
          #                      -------------------------------------
          #               locus1           >>>>>>>>>>
          # or ===============================================================
          #               locus2                 <<<<<<<<<
          #                      -------------------------------------
          #               locus1                     <<<<<<<<<

            ###    3           if locus1 5' is inside locus2, but locus1 3' extends locus2 3'

          #               locus2                 >>>>>>>>>
          #                      -------------------------------------
          #               locus1                     >>>>>>>>>>
          # or ===============================================================
          #               locus2                 <<<<<<<<<
          #                      -------------------------------------
          #               locus1              <<<<<<<<<

            ###    4           if locus2 is fully included in locus1

          #               locus2                 >>>>>>>>>
          #                      -------------------------------------
          #               locus1             >>>>>>>>>>>>>>>>>>
          # or ===============================================================
          #               locus2                 <<<<<<<<<
          #                      -------------------------------------
          #               locus1             <<<<<<<<<<<<<<<<<<

            ###    5           if locus1 5' is strictly inside locus2, but locus1 and locus2 share the same 3' end

          #               locus2                 >>>>>>>>>
          #                      -------------------------------------
          #               locus1                     >>>>>
          # or ===============================================================
          #               locus2                 <<<<<<<<<
          #                      -------------------------------------
          #               locus1                 <<<<<

            ###    6           if locus1 3' is strictly inside locus2, but locus1 and locus2 share the same 5' end

          #               locus2                 >>>>>>>>>
          #                      -------------------------------------
          #               locus1                 >>>>>
          # or ===============================================================
          #               locus2                 <<<<<<<<<
          #                      -------------------------------------
          #               locus1                     <<<<<

            ###    7           if locus1 and locus2 share the exact same coordinates

          #               locus2                 >>>>>>>>>
          #                      -------------------------------------
          #               locus1                 >>>>>>>>>
          # or ===============================================================
          #               locus2                 <<<<<<<<<
          #                      -------------------------------------
          #               locus1                 <<<<<<<<<

            ###    8           if locus1's 5'end extends locus2's, but locus1 and locus2 share the same 3' end

          #               locus2                 >>>>>>>>>
          #                      -------------------------------------
          #               locus1            >>>>>>>>>>>>>>
          # or ===============================================================
          #               locus2                 <<<<<<<<<
          #                      -------------------------------------
          #               locus1                 <<<<<<<<<<<<<<<

#             9 if locus1's 3'end extends locus2's, but locus1 and locus2 share the same 5' end
#               locus2                 >>>>>>>>>
#                      -------------------------------------
#               locus1                 >>>>>>>>>>>>>
# or ===============================================================
#               locus2                 <<<<<<<<<
#                      -------------------------------------
#               locus1             <<<<<<<<<<<<<

            if ( $locus1chr eq $locus2chr && $locus1strand eq $locus2strand )
            {
                # possible values -2,-1,1,2,3,4
                my $fiveprimeDist;
                my $threeprimeDist;
                if ( $locus1strand eq '+' ) {
                    $fiveprimeDist  = $startlocus1 - $startlocus2;
                    $threeprimeDist = $endlocus1 - $endlocus2;
                }
                else {
                    $fiveprimeDist  = $endlocus2 - $endlocus1;
                    $threeprimeDist = $startlocus2 - $startlocus1;
                }

                # enumerating cases one by one, from left to right

                if ( $startlocus1 > $endlocus2 ) {    #no overlap
                                                      #possible values -2,-1
                    if ( $locus1strand eq '+' ) {
                        return ( -1, $fiveprimeDist, $threeprimeDist );
                    }
                    else {
                        return ( -2, $fiveprimeDist, $threeprimeDist );
                    }
                }

                else {    # i.e. ($startlocus1 <= $endlocus2)
                    if ( $endlocus1 < $startlocus2 )
                    { # i.e. ($startlocus1 <= $endlocus2) && ($endlocus1 < $startlocus2) # no overlap
                        if ( $locus1strand eq '+' ) {
                            return ( -2, $fiveprimeDist, $threeprimeDist );
                        }
                        else {
                            return ( -1, $fiveprimeDist, $threeprimeDist );
                        }
                    }
                    elsif ( $startlocus1 < $startlocus2 ) {
                        if ( $endlocus1 >= $startlocus2 ) {
                            if ( $endlocus1 < $endlocus2 ) {
                                if ( $locus1strand eq '+' ) {
                                    return ( 2, $fiveprimeDist,
                                        $threeprimeDist );
                                }
                                else {
                                    return ( 3, $fiveprimeDist,
                                        $threeprimeDist );
                                }
                            }
                            elsif ( $endlocus1 == $endlocus2 ) {
                                if ( $locus1strand eq '+' ) {
                                    return ( 8, $fiveprimeDist,
                                        $threeprimeDist );
                                }
                                else {
                                    return ( 9, $fiveprimeDist,
                                        $threeprimeDist );
                                }
                            }
                            elsif ( $endlocus1 > $endlocus2 ) {
                                return ( 4, $fiveprimeDist, $threeprimeDist );
                            }
                            else {
                                die;
                            }
                        }
                    }
                    elsif ( $startlocus1 == $startlocus2 ) {
                        if ( $endlocus1 < $endlocus2 ) {
                            if ( $locus1strand eq '+' ) {
                                return ( 6, $fiveprimeDist, $threeprimeDist );
                            }
                            else {
                                return ( 5, $fiveprimeDist, $threeprimeDist );
                            }
                        }
                        elsif ( $endlocus1 == $endlocus2 ) {
                            return ( 7, $fiveprimeDist, $threeprimeDist );
                        }
                        elsif ( $endlocus1 > $endlocus2 ) {
                            if ( $locus1strand eq '+' ) {
                                return ( 9, $fiveprimeDist, $threeprimeDist );
                            }
                            else {
                                return ( 8, $fiveprimeDist, $threeprimeDist );
                            }
                        }
                        else {
                            die;
                        }
                    }
                    elsif ( $startlocus1 > $startlocus2 ) {
                        if ( $endlocus1 < $endlocus2 ) {
                            return ( 1, $fiveprimeDist, $threeprimeDist );
                        }
                        elsif ( $endlocus1 == $endlocus2 ) {
                            if ( $locus1strand eq '+' ) {
                                return ( 5, $fiveprimeDist, $threeprimeDist );
                            }
                            else {
                                return ( 6, $fiveprimeDist, $threeprimeDist );
                            }
                        }
                        elsif ( $endlocus1 > $endlocus2 ) {
                            if ( $startlocus1 <= $endlocus2 ) {
                                if ( $locus1strand eq '+' ) {
                                    return ( 3, $fiveprimeDist,
                                        $threeprimeDist );
                                }
                                else {
                                    return ( 2, $fiveprimeDist,
                                        $threeprimeDist );
                                }
                            }
                            else {
                                die;
                            }
                        }
                        else {
                            die;
                        }
                    }
                    else {
                        die;
                    }
                }
            }
            else {
                return ( 0, undef, undef );
            }
        }

        1;
