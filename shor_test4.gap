#############################################################################
# MINSPM -- Dictionary-Based QFT Period-Finding Simulator
#############################################################################

MINSPM_Gcd := function(a, b)
    local tmp;

    a := AbsInt(a);
    b := AbsInt(b);

    while b <> 0 do
        tmp := a mod b;
        a := b;
        b := tmp;
    od;

    return a;
end;


MINSPM_PowerMod := function(a, e, n)
    local result, base, exp;

    result := 1 mod n;
    base := a mod n;
    exp := e;

    while exp > 0 do
        if (exp mod 2) = 1 then
            result := (result * base) mod n;
        fi;

        base := (base * base) mod n;
        exp := QuoInt(exp, 2);
    od;

    return result;
end;


#############################################################################
# Build the actual modular-exponentiation dictionary
#
# Each entry is:
#
#     [ f(x), multiplicity ]
#
# where
#
#     f(x) = a^x mod N
#
#############################################################################

MINSPMBuildModExpDictionary := function(Q, N, a)
    local x, value, key, rec_dict, dict, names;

    rec_dict := rec();

    for x in [0 .. Q - 1] do
        value := MINSPM_PowerMod(a, x, N);
        key := String(value);

        if IsBound(rec_dict.(key)) then
            rec_dict.(key) := rec_dict.(key) + 1;
        else
            rec_dict.(key) := 1;
        fi;
    od;

    dict := [];
    names := RecNames(rec_dict);

    for key in names do
        Add(dict, [Int(key), rec_dict.(key)]);
    od;

    Sort(dict, function(x, y)
        return x[1] < y[1];
    end);

    return dict;
end;


#############################################################################
# Find the x-values associated with a particular measured modular value.
#
# This represents the post-measurement first-register state.
#############################################################################

MINSPMBuildCoset := function(Q, N, a, measured_value)
    local x, coset;

    coset := [];

    for x in [0 .. Q - 1] do
        if MINSPM_PowerMod(a, x, N) = measured_value then
            Add(coset, x);
        fi;
    od;

    return coset;
end;


#############################################################################
# Complex arithmetic
#############################################################################

MINSPM_ComplexAdd := function(a, b)
    return [a[1] + b[1], a[2] + b[2]];
end;


MINSPM_ComplexScale := function(c, a)
    return [c * a[1], c * a[2]];
end;


MINSPM_ComplexNormSquared := function(a)
    return a[1] * a[1] + a[2] * a[2];
end;


MINSPM_PI := 3.1415926535897932384626433832795;


#############################################################################
# QFT amplitude for a particular output frequency y
#
# alpha_y =
#
#   1/sqrt(Q L) * Sum_x exp(2*pi*i*x*y/Q)
#
# over the measured coset.
#############################################################################

MINSPM_QFTAmplitude := function(coset, Q, y)
    local result, x, angle, L;

    L := Length(coset);

    if L = 0 then
        return [0.0, 0.0];
    fi;

    result := [0.0, 0.0];

    for x in coset do
        angle := 2.0 * MINSPM_PI * x * y / Q;

        result[1] := result[1] + Cos(angle);
        result[2] := result[2] + Sin(angle);
    od;

    return MINSPM_ComplexScale(
        1.0 / Sqrt(Q * L),
        result
    );
end;


MINSPM_QFTProbability := function(coset, Q, y)
    local z;

    z := MINSPM_QFTAmplitude(coset, Q, y);

    return MINSPM_ComplexNormSquared(z);
end;


#############################################################################
# Compute the complete QFT probability distribution.
#
# This is O(Q * L), so it is a classical simulation rather than a
# quantum implementation.
#############################################################################

MINSPM_QFTSpectrum := function(coset, Q)
    local spectrum, y;

    spectrum := [];

    for y in [0 .. Q - 1] do
        Add(
            spectrum,
            MINSPM_QFTProbability(coset, Q, y)
        );
    od;

    return spectrum;
end;


#############################################################################
# Score a candidate period r.
#
# We do not merely look at y=1.
#
# For candidate r, expected peaks satisfy approximately
#
#       y/Q ~= s/r.
#
# We therefore collect probability around the nearest predicted
# frequencies.
#############################################################################

MINSPM_PeriodScore := function(coset, Q, r)
    local score, s, y, width, center, yy, y0, y1;

    if r < 1 then
        return 0.0;
    fi;

    width := 2;
    score := 0.0;

    for s in [0 .. r - 1] do

        center := s * Q / r;
        y := Int(center + 0.5);

        y0 := Maximum(0, y - width);
        y1 := Minimum(Q - 1, y + width);

        for yy in [y0 .. y1] do
            score := score +
                MINSPM_QFTProbability(coset, Q, yy);
        od;
    od;

    return score;
end;


#############################################################################
# Verify an actual multiplicative order candidate.
#############################################################################

MINSPM_VerifyPeriod := function(N, a, r)
    if r <= 0 then
        return false;
    fi;

    if MINSPM_Gcd(a, N) <> 1 then
        return false;
    fi;

    return MINSPM_PowerMod(a, r, N) = 1;
end;


#############################################################################
# Search candidate periods.
#############################################################################

MINSPMOptimizePeriod := function(N, base_a, estimated_start, search_radius)
    local Q, dict, measured_value, coset;
    local best_r, best_score, r, score;
    local candidates, item;

    Q := 2^16;

    Print("\n[MINSPM] Building modular-exponentiation dictionary...\n");

    dict := MINSPMBuildModExpDictionary(Q, N, base_a);

    Print(
        "[MINSPM] Distinct modular values: ",
        Length(dict), "\n"
    );

    #
    # In a real period-finding experiment the second register is measured.
    # For simulation we select one value from the dictionary.
    #
    measured_value := dict[1][1];

    coset := MINSPMBuildCoset(
        Q,
        N,
        base_a,
        measured_value
    );

    Print(
        "[MINSPM] Measured value: ",
        measured_value,
        " | coset size: ",
        Length(coset),
        "\n"
    );

    best_r := 0;
    best_score := -1.0;

    candidates := [
        Maximum(1, estimated_start - search_radius) ..
        estimated_start + search_radius
    ];

    for r in candidates do

        score := MINSPM_PeriodScore(
            coset,
            Q,
            r
        );

        if score > best_score then
            best_score := score;
            best_r := r;
        fi;
    od;

    Print(
        "[MINSPM] Best spectral candidate: r = ",
        best_r,
        " | score = ",
        best_score,
        "\n"
    );

    if MINSPM_VerifyPeriod(N, base_a, best_r) then
        Print(
            "[MINSPM] Candidate verified: ",
            base_a, "^", best_r,
            " = 1 (mod ", N, ")\n"
        );
    else
        Print(
            "[MINSPM] Candidate is NOT the multiplicative order.\n"
        );
    fi;

    return best_r;
end;


#############################################################################
# Example
#############################################################################

Target_N := 10403;
Target_a := 2;
Estimated_Start := 45;

Found_Period :=
    MINSPMOptimizePeriod(
        Target_N,
        Target_a,
        Estimated_Start,
        20
    );

Print(
    "\nSpectral Candidate Period: ",
    Found_Period,
    "\n"
);

