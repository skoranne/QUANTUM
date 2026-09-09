#############################################################################
#
# MINSPM + AQFT PERIOD FINDING / FACTORIZATION
#
# Self-contained GAP 4 program.
#
# The MINSPM dictionary stores sparse Fourier exponents:
#
#       A(k) = Sum_j c_j * exp(2*pi*i*e_j/Q)
#
# where Q = 2^t.
#
# IMPORTANT:
# We do NOT identify e with -e.  Those are distinct phases.
#
# The AQFT simulator uses the standard order-finding conditional state
#
#       |psi> = 1/sqrt(L) Sum_j |x0 + j*r>
#
# followed by the inverse QFT:
#
#       A(k) = 1/sqrt(Q*L)
#              Sum_j exp(2*pi*i*k*(x0+j*r)/Q).
#
# Only the sparse phase dictionary is evaluated, so we do not construct
# a Q x Q Fourier matrix.
#
#############################################################################


#############################################################################
# Basic integer utilities
#############################################################################

MINSPM_Gcd := function(a,b)
    while b <> 0 do
        a := a mod b;
        if a = 0 then
            return AbsInt(b);
        fi;
        # Exchange using a temporary.
        # GAP has no tuple assignment.
        tmp := a;
        a := b;
        b := tmp;
    od;
    return AbsInt(a);
end;


MINSPM_GcdExtended := function(a,b)
    local r0,r1,s0,s1,t0,t1,q,tmp_r,tmp_s,tmp_t;

    r0 := a;
    r1 := b;

    s0 := 1;
    s1 := 0;

    t0 := 0;
    t1 := 1;

    while r1 <> 0 do
        q := QuoInt(r0,r1);

        tmp_r := r1;
        r1 := r0 - q*r1;
        r0 := tmp_r;

        tmp_s := s1;
        s1 := s0 - q*s1;
        s0 := tmp_s;

        tmp_t := t1;
        t1 := t0 - q*t1;
        t0 := tmp_t;
    od;

    return rec(gcd := AbsInt(r0), x := s0, y := t0);
end;


#############################################################################
# Modular arithmetic
#############################################################################

MINSPM_PowerMod := function(a,e,n)
    local result,base,exp;

    result := 1 mod n;
    base := a mod n;
    exp := e;

    while exp > 0 do
        if (exp mod 2) = 1 then
            result := (result * base) mod n;
        fi;

        base := (base * base) mod n;
        exp := QuoInt(exp,2);
    od;

    return result;
end;


#############################################################################
# Classical order routine.
#
# This routine is ONLY used by the simulator to prepare the synthetic
# quantum/AQFT measurement distribution.  It is NOT used to recover the
# period from the measurement.
#############################################################################

MINSPMClassicalOrder := function(a,n)
    local x,r;

    if MINSPM_Gcd(a,n) <> 1 then
        return fail;
    fi;

    x := 1;
    r := 0;

    repeat
        x := (x*a) mod n;
        r := r + 1;

        if r > n then
            return fail;
        fi;
    until x = 1;

    return r;
end;


#############################################################################
# MINSPM sparse phase dictionary
#############################################################################

# NOTE:
# There is deliberately NO canonical +/- exponent folding here.
#
# E(Q)^e and E(Q)^(-e) are different complex phases.
#
# The dictionary is represented as:
#
#       [ [ exponent1, coefficient1 ],
#         [ exponent2, coefficient2 ],
#         ... ]
#
# with distinct exponents.

MINSPMBuildDictionary := function(exponents,Q)
    local work,dict,i,e,found,j;

    work := [];

    for e in exponents do
        Add(work,[e mod Q,1]);
    od;

    if Length(work) = 0 then
        return [];
    fi;

    Sort(work,
        function(x,y)
            return x[1] < y[1];
        end);

    dict := [];

    for i in [1..Length(work)] do
        e := work[i][1];

        if Length(dict) = 0 then
            Add(dict,[e,work[i][2]]);
        else
            j := Length(dict);

            if dict[j][1] = e then
                dict[j][2] := dict[j][2] + work[i][2];
            else
                Add(dict,[e,work[i][2]]);
            fi;
        fi;
    od;

    return dict;
end;


#############################################################################
# Floating-point complex pair arithmetic
#
# A complex number is represented by [real,imaginary].
#############################################################################

MINSPM_ComplexAdd := function(a,b)
    return [a[1]+b[1],a[2]+b[2]];
end;


MINSPM_ComplexScale := function(c,a)
    return [c*a[1],c*a[2]];
end;


MINSPM_ComplexNormSquared := function(a)
    return a[1]*a[1] + a[2]*a[2];
end;


#############################################################################
# Dyadic phase
#############################################################################

MINSPM_PI := 3.1415926535897932384626433832795;


MINSPMPhaseFloat := function(e,Q)
    local angle;

    angle := 2.0 * MINSPM_PI * (e mod Q) / Q;

    return [Cos(angle),Sin(angle)];
end;


#############################################################################
# Evaluate sparse MINSPM dictionary
#############################################################################

MINSPMEvaluateDictionary := function(dict,Q)
    local result,item,phase,term;

    result := [0.0,0.0];

    for item in dict do
        phase := MINSPMPhaseFloat(item[1],Q);
        term := MINSPM_ComplexScale(item[2],phase);
        result := MINSPM_ComplexAdd(result,term);
    od;

    return result;
end;


#############################################################################
# Generate the actual order-finding Fourier exponents.
#
# Conditional state:
#
#     x = x0 + j*r
#
# for all x<Q.
#
# After QFT:
#
#     phase exponent = k*x mod Q.
#############################################################################

MINSPMGenerateAQFTExponents := function(Q,r,x0,k)
    local exponents,x,j,L;

    exponents := [];

    if x0 >= Q then
        x0 := x0 mod r;
    fi;

    j := 0;
    x := x0;

    while x < Q do
        Add(exponents,(k*x) mod Q);

        j := j + 1;
        x := x0 + j*r;
    od;

    return exponents;
end;


#############################################################################
# AQFT spectral amplitude
#############################################################################

MINSPMAQFTAmplitude := function(t,r,x0,k)
    local Q,exponents,dict,z,L;

    Q := 2^t;

    exponents := MINSPMGenerateAQFTExponents(Q,r,x0,k);

    L := Length(exponents);

    if L = 0 then
        return [0.0,0.0];
    fi;

    dict := MINSPMBuildDictionary(exponents,Q);

    z := MINSPMEvaluateDictionary(dict,Q);

    # Normalization:
    #
    # inverse QFT contributes 1/sqrt(Q)
    # conditional state contributes 1/sqrt(L)
    #
    return MINSPM_ComplexScale(
        1.0 / Sqrt(Q*L),
        z
    );
end;


#############################################################################
# AQFT spectral intensity
#############################################################################

MINSPMAQFTIntensity := function(t,r,x0,k)
    local amp;

    amp := MINSPMAQFTAmplitude(t,r,x0,k);

    return MINSPM_ComplexNormSquared(amp);
end;


#############################################################################
# MINSPM spectral profile around a theoretical AQFT peak
#############################################################################

MINSPMProfileAroundPeak := function(t,r,x0,m,window)
    local Q,k0,k,profile;

    Q := 2^t;

    # k0 ~= m*Q/r
    k0 := QuoInt(m*Q + QuoInt(r,2),r);

    profile := [];

    for k in [Maximum(0,k0-window)..Minimum(Q-1,k0+window)] do
        Add(profile,
            [k,MINSPMAQFTIntensity(t,r,x0,k)]);
    od;

    return profile;
end;


#############################################################################
# Continued-fraction convergents
#############################################################################

MINSPMContinuedFractionConvergents := function(num,den)
    local n,d,a,p0,p1,q0,q1,p2,q2,tmp,ans;

    ans := [];

    if den = 0 then
        return ans;
    fi;

    # Reduce the rational number first.
    tmp := MINSPM_Gcd(AbsInt(num),AbsInt(den));

    if tmp <> 0 then
        num := QuoInt(num,tmp);
        den := QuoInt(den,tmp);
    fi;

    n := num;
    d := den;

    p0 := 0;
    p1 := 1;

    q0 := 1;
    q1 := 0;

    while d <> 0 do
        a := QuoInt(n,d);

        p2 := a*p1 + p0;
        q2 := a*q1 + q0;

        Add(ans,[p2,q2]);

        p0 := p1;
        p1 := p2;

        q0 := q1;
        q1 := q2;

        tmp := n - a*d;
        n := d;
        d := tmp;
    od;

    return ans;
end;


#############################################################################
# Recover candidate order from an AQFT measurement k/Q.
#############################################################################

MINSPMRecoverOrder := function(a,n,k,Q)
    local conv,cand,p,q,seen;

    conv := MINSPMContinuedFractionConvergents(k,Q);

    seen := [];

    for cand in conv do
        p := cand[1];
        q := cand[2];

        if q > 1 and q <= n then

            if not q in seen then
                Add(seen,q);

                if MINSPM_PowerMod(a,q,n) = 1 then
                    return q;
                fi;
            fi;
        fi;
    od;

    return fail;
end;


#############################################################################
# Find a non-trivial Shor factor from an even order.
#############################################################################

MINSPMSplitFactor := function(a,n,r)
    local x,g;

    if r = fail then
        return fail;
    fi;

    if r mod 2 <> 0 then
        return fail;
    fi;

    x := MINSPM_PowerMod(a,QuoInt(r,2),n);

    if x = 1 or x = n-1 then
        return fail;
    fi;

    g := MINSPM_Gcd(x-1,n);

    if g > 1 and g < n then
        return [g,QuoInt(n,g)];
    fi;

    return fail;
end;


#############################################################################
# Simulated AQFT measurement.
#
# In an ideal order-finding circuit, peaks occur near
#
#       k/Q = m/r.
#
# We choose a coprime m and then select the nearest dyadic output k.
#
# The resulting k is subsequently fed through the SAME continued-fraction
# recovery procedure that would process a measured quantum result.
#
# The MINSPM dictionary is evaluated at the selected k so that the program
# actually computes the corresponding AQFT spectral amplitude.
#############################################################################

MINSPMSampleAQFT := function(t,r,x0)
    local Q,m,k,intensity;

    Q := 2^t;

    # Avoid a reducible m/r.
    repeat
        m := Random(1,r-1);
    until MINSPM_Gcd(m,r) = 1;

    # Nearest integer to m*Q/r.
    k := QuoInt(m*Q + QuoInt(r,2),r);

    if k >= Q then
        k := Q-1;
    fi;

    intensity := MINSPMAQFTIntensity(t,r,x0,k);

    return rec(
        k := k,
        Q := Q,
        m := m,
        intensity := intensity,
        ratio := [k,Q]
    );
end;


#############################################################################
# Complete MINSPM/AQFT factoring driver
#############################################################################

FactorWithMINSPMAQFT := function(N,base_a)
    local t_gates,Q,oracle_r,x0,sample,recovered_r,factors,
          attempt,max_attempts,start_time;

    # For order r <= N, the usual continued-fraction condition is
    #
    #       |k/Q - m/r| < 1/(2*r^2).
    #
    # Since Q is dyadic, choose Q > 2*N^2.
    #
    # 2^27 = 134217728 > 2*(10403)^2.
    #
    t_gates := 27;
    Q := 2^t_gates;

    Print("\n");
    Print("============================================================\n");
    Print(" MINSPM + AQFT ORDER-FINDING FACTORIZATION\n");
    Print("============================================================\n");
    Print("N              = ",N,"\n");
    Print("base a         = ",base_a,"\n");
    Print("AQFT qubits    = ",t_gates,"\n");
    Print("AQFT dimension = ",Q,"\n");
    Print("============================================================\n");

    if MINSPM_Gcd(base_a,N) <> 1 then
        Print("gcd(a,N) is already non-trivial.\n");
        return [
            MINSPM_Gcd(base_a,N),
            QuoInt(N,MINSPM_Gcd(base_a,N))
        ];
    fi;

    start_time := Runtime();

    #########################################################################
    # Simulation-only step:
    #
    # We need the actual order to prepare the synthetic quantum state.
    #
    # This value is NOT passed to the recovery routine.
    #########################################################################

    oracle_r := MINSPMClassicalOrder(base_a,N);

    if oracle_r = fail then
        Print("Unable to determine simulator order.\n");
        return fail;
    fi;

    Print("Simulator order = ",oracle_r,"\n");
    Print("(Used only to generate synthetic AQFT measurements.)\n");
    Print("------------------------------------------------------------\n");

    max_attempts := 20;

    for attempt in [1..max_attempts] do

        # The measured second register selects an offset x0.
        x0 := Random(0,oracle_r-1);

        sample := MINSPMSampleAQFT(
            t_gates,
            oracle_r,
            x0
        );

        Print("AQFT shot ",attempt,
              ": k = ",sample.k,
              ", Q = ",sample.Q,
              ", k/Q = ",
              sample.k,"/",sample.Q,
              ", peak intensity = ",
              sample.intensity,"\n");

        #####################################################################
        # Classical post-processing of the quantum measurement.
        #####################################################################

        recovered_r := MINSPMRecoverOrder(
            base_a,
            N,
            sample.k,
            sample.Q
        );

        if recovered_r <> fail then

            Print("  -> Continued fractions recovered r = ",
                  recovered_r,"\n");

            if MINSPM_PowerMod(base_a,recovered_r,N) = 1 then

                factors := MINSPMSplitFactor(
                    base_a,
                    N,
                    recovered_r
                );

                if factors <> fail then
                    Print("------------------------------------------------------------\n");
                    Print("SUCCESS\n");
                    Print("Recovered order = ",recovered_r,"\n");
                    Print("Factor 1        = ",factors[1],"\n");
                    Print("Factor 2        = ",factors[2],"\n");
                    Print("Product         = ",
                          factors[1]*factors[2],"\n");
                    Print("Runtime         = ",
                          (Runtime()-start_time)/1000.0,
                          " seconds\n");
                    Print("============================================================\n");

                    return rec(
                        factors := factors,
                        order := recovered_r,
                        measurement := sample,
                        qubits := t_gates
                    );
                fi;

                Print("  -> Order is valid but Shor split was trivial.\n");

            fi;
        else
            Print("  -> Continued fractions did not recover the order.\n");
        fi;
    od;

    Print("------------------------------------------------------------\n");
    Print("No non-trivial factor obtained in ",max_attempts," AQFT shots.\n");
    Print("Try increasing the number of shots or t_gates.\n");
    Print("============================================================\n");

    return fail;
end;


#############################################################################
# Optional diagnostic: print the strongest points around an AQFT peak.
#############################################################################

MINSPMShowPeakProfile := function(N,base_a,t,m,x0)
    local r,Q,profile,item;

    r := MINSPMClassicalOrder(base_a,N);
    Q := 2^t;

    Print("\nMINSPM/AQFT peak profile\n");
    Print("N = ",N,", a = ",base_a,", r = ",r,", Q = ",Q,"\n");
    Print("m = ",m,", x0 = ",x0,"\n\n");

    profile := MINSPMProfileAroundPeak(
        t,
        r,
        x0,
        m,
        5
    );

    for item in profile do
        Print("k = ",item[1],
              "    intensity = ",
              item[2],"\n");
    od;

    return profile;
end;


#############################################################################
# Example
#############################################################################

Target_N := 10403;
Target_a := 2;

result := FactorWithMINSPMAQFT(
    Target_N,
    Target_a
);

Print("\nReturned object:\n");
Print(result,"\n");

