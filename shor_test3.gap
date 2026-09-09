#############################################################################
# MINSPM + AQFT PERIOD FINDING / FACTORIZATION (OPTIMIZED)
# Self-contained GAP 4 program.
#
# The MINSPM dictionary stores sparse Fourier exponents via hash-based
# frequency accumulation, avoiding massive array allocation and sorting overhead.
#############################################################################

MINSPM_Gcd := function(a,b)
    local tmp;
    while b <> 0 do
        a := a mod b;
        if a = 0 then
            return AbsInt(b);
        fi;
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

# Optimized O(L) Dictionary Builder using Record Hash Accumulation
MINSPMBuildDictionaryOptimized := function(Q, r, x0, k)
    local A, B, j, x, e, rec_dict, names, key, dict;
    
    if x0 >= Q then
        x0 := x0 mod r;
    fi;

    A := (k * x0) mod Q;
    B := (k * r) mod Q;
    
    rec_dict := rec();
    
    j := 0;
    x := x0;
    while x < Q do
        e := (A + j * B) mod Q;
        key := String(e);
        if IsBound(rec_dict.(key)) then
            rec_dict.(key) := rec_dict.(key) + 1;
        else
            rec_dict.(key) := 1;
        fi;
        
        j := j + 1;
        x := x0 + j * r;
    od;
    
    dict := [];
    names := RecNames(rec_dict);
    for key in names do
        e := Int(key);
        Add(dict, [e, rec_dict.(key)]);
    od;
    
    Sort(dict, function(x, y) return x[1] < y[1]; end);
    return dict;
end;

MINSPM_ComplexAdd := function(a,b)
    return [a[1]+b[1],a[2]+b[2]];
end;

MINSPM_ComplexScale := function(c,a)
    return [c*a[1],c*a[2]];
end;

MINSPM_ComplexNormSquared := function(a)
    return a[1]*a[1] + a[2]*a[2];
end;

MINSPM_PI := 3.1415926535897932384626433832795;

MINSPMPhaseFloat := function(e,Q)
    local angle;
    angle := 2.0 * MINSPM_PI * (e mod Q) / Q;
    return [Cos(angle),Sin(angle)];
end;

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

MINSPMAQFTAmplitude := function(t,r,x0,k)
    local Q,dict,z,L,item;
    Q := 2^t;
    
    dict := MINSPMBuildDictionaryOptimized(Q, r, x0, k);

    L := 0;
    for item in dict do
        L := L + item[2];
    od;

    if L = 0 then
        return [0.0,0.0];
    fi;

    z := MINSPMEvaluateDictionary(dict,Q);

    return MINSPM_ComplexScale(
        1.0 / Sqrt(Q*L),
        z
    );
end;

MINSPMAQFTIntensity := function(t,r,x0,k)
    local amp;
    amp := MINSPMAQFTAmplitude(t,r,x0,k);
    return MINSPM_ComplexNormSquared(amp);
end;

MINSPMProfileAroundPeak := function(t,r,x0,m,window)
    local Q,k0,k,profile;
    Q := 2^t;
    k0 := QuoInt(m*Q + QuoInt(r,2),r);
    profile := [];
    for k in [Maximum(0,k0-window)..Minimum(Q-1,k0+window)] do
        Add(profile, [k, MINSPMAQFTIntensity(t,r,x0,k)]);
    od;
    return profile;
end;

MINSPMContinuedFractionConvergents := function(num,den)
    local n,d,a,p0,p1,q0,q1,p2,q2,tmp,ans;
    ans := [];
    if den = 0 then
        return ans;
    fi;
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

MINSPMSplitFactor := function(a,n,r)
    local x,g;
    if r = fail or r mod 2 <> 0 then
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

MINSPMSampleAQFT := function(t,r,x0)
    local Q,m,k,intensity;
    Q := 2^t;
    repeat
        m := Random(1,r-1);
    until MINSPM_Gcd(m,r) = 1;
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

FactorWithMINSPMAQFT := function(N,base_a)
    local t_gates,Q,oracle_r,x0,sample,recovered_r,factors,
          attempt,max_attempts,start_time;

    t_gates := 27;
    Q := 2^t_gates;

    Print("\n");
    Print("============================================================\n");
    Print(" MINSPM + AQFT ORDER-FINDING FACTORIZATION (OPTIMIZED)\n");
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
    oracle_r := MINSPMClassicalOrder(base_a,N);

    if oracle_r = fail then
        Print("Unable to determine simulator order.\n");
        return fail;
    fi;

    Print("Simulator order = ",oracle_r,"\n");
    Print("------------------------------------------------------------\n");

    max_attempts := 20;

    for attempt in [1..max_attempts] do
        x0 := Random(0,oracle_r-1);
        sample := MINSPMSampleAQFT(t_gates, oracle_r, x0);

        Print("AQFT shot ",attempt,
              ": k = ",sample.k,
              ", Q = ",sample.Q,
              ", k/Q = ",
              sample.k,"/",sample.Q,
              ", peak intensity = ",
              sample.intensity,"\n");

        recovered_r := MINSPMRecoverOrder(
            base_a,
            N,
            sample.k,
            sample.Q
        );

        if recovered_r <> fail then
            Print("  -> Continued fractions recovered r = ", recovered_r,"\n");

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
                    Print("Product         = ",factors[1]*factors[2],"\n");
                    Print("Runtime         = ",(Runtime()-start_time)/1000.0," seconds\n");
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
    Print("============================================================\n");

    return fail;
end;

Target_N := 10403;
Target_a := 2;

result := FactorWithMINSPMAQFT(
    Target_N,
    Target_a
);

Print("\nReturned object:\n");
Print(result,"\n");
