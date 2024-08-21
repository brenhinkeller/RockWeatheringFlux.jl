"""
```julia 
screen_outliers!(bulk, cats; [warn])
```

Remove geochemical outliers from `bulk`. Some outliers are removed by rock class.
Specify rock classes in `cats` (see `match_rocktype`). Optionally set `warn=true` to print 
a warning to the terminal if an expected element or element oxide is not present in `bulk`.

"""
function screen_outliers!(bulk::NamedTuple, cats::NamedTuple; warn::Bool=false)
    # Corrections are based on outlier screening for EarthChem (EC values, done by 
    # Brenhin) and screening specifically for the Gard dataset (done by Brenhin and Blair).
    # Note both screenings were done for igneous rocks only. 
    # 
    # All lower bounds have been adjusted to 0; all upper bounds have been re-checked
    # against a histogram of all rock class and the upper 99th percentile for that 
    # element (for each rock class). All values have also been converted to wt.%: 
    # INCLUDING trace elements.
    #
    # The source of values is noted as EC (EarthChem) Gard (Gard) or histogram / percentile
    # (99th [class])

    # Cl, F, S are bimodal by four orders of magnitude. We converted everything 
    # ppm -> wt.% (./10_000). Assume anything under 0.0003% (3 ppm) was in wt.% and 
    # needs un-converted. From Gard and EarthChem screening
    bulk.Cl[bulk.Cl .< 0.0003] .= bulk.Cl[bulk.Cl .< 0.0003] .* 10_000;
    bulk.F[bulk.F .< 0.0003] .= bulk.F[bulk.F .< 0.0003] .* 10_000;
    bulk.S[bulk.S .< 0.0003] .= bulk.S[bulk.S .< 0.0003] .* 10_000;

    # Some data were actually PPB, so they're too large and need an additional division
    # From Gard screening
    bulk.Au[bulk.Au .> 5e-5] .= bulk.Au[bulk.Au .> 5e-5] ./ 1000
    bulk.Bi[bulk.Bi .> 0.001] .= bulk.Bi[bulk.Bi .> 0.001] ./ 1000
    bulk.Pd[bulk.Pd .> 1e-5] .= bulk.Pd[bulk.Pd .> 1e-5] ./ 1000
    bulk.Hg[bulk.Hg .> 0.001] .= bulk.Hg[bulk.Hg .> 0.001] ./ 1000
    bulk.Pt[bulk.Pt .> 1e-5] .= bulk.Pt[bulk.Pt .> 1e-5] ./ 1000
    bulk.W[bulk.W .> 0.001] .= bulk.W[bulk.W .> 0.001] ./ 1000;

    # Same hat but I found these myself and want validation about my cutoffs 
    bulk.Ag[bulk.Ag .> 0.005] .= bulk.Ag[bulk.Ag .> 0.005] ./ 1000;     # Feels ok
    bulk.As[bulk.As .> 0.005] .= bulk.As[bulk.As .> 0.005] ./ 1000;     # Feels good

    # Weird that there's so many that are exactly the same evenly spaced value !!
    for i in 0.0001:0.0001:0.001
        bulk.W[bulk.W .== i] .= NaN
    end

    # These bad boys could have been PPT, PPM, or PPB: histograms are bunched in the 
    # smallest bins with really large 99th percentiles
    bulk.Re[bulk.Re .> 0.1] .= bulk.Re[bulk.Re .> 0.1] ./ 1e9
    bulk.Re[bulk.Re .> 0.001] .= bulk.Re[bulk.Re .> 0.001] ./ 1e6
    bulk.Re[bulk.Re .> 1e-6] .= bulk.Re[bulk.Re .> 1e-6] ./ 1e3

    bulk.Os[bulk.Os .> 0.1] .= bulk.Os[bulk.Os .> 0.1] ./ 1e9
    bulk.Os[bulk.Os .> 0.001] .= bulk.Os[bulk.Os .> 0.001] ./ 1e6
    bulk.Os[bulk.Os .> 1e-6] .= bulk.Os[bulk.Os .> 1e-6] ./ 1e3

    bulk.Ir[bulk.Ir .> 0.1] .= bulk.Ir[bulk.Ir .> 0.1] ./ 1e9
    bulk.Ir[bulk.Ir .> 0.001] .= bulk.Ir[bulk.Ir .> 0.001] ./ 1e6
    bulk.Ir[bulk.Ir .> 1e-6] .= bulk.Ir[bulk.Ir .> 1e-6] ./ 1e3

    # Other elements that may have mixed units: Cd, Cl, Co, Cs, Cu, I, In, Mo, Pb, 
    # Sb, Se, Sn, Ta, Te, Tl, U, Zn, (P₂O₅ ?)
    
    # Reasonable min / max values for each element. Min is always 0, for fun
    bounds = (
        # Major element oxides 
        SiO2=(0,100), Al2O3=(0,30),
        # FeOT=(0,30),
        TiO2=(0,10),
        MgO=(0,50),     # Forsterite (Mg₂SiO₄) = 27% MgO; 99th met = 42%. Gard has 50%
        CaO=(0,56),     # Calcite (CaCO₃) = 56% CaO. 99th sed = 55%
        Na2O=(0,20),    # EC. 30, reduced from histogram to match Gard. 99th = 7.5% 
        K2O=(0,30),

        # Minor elements and element oxides
        Ag=(0,0.001),   # EC 0.3, changed to match Gard; 99th met = 0.04
        As=(0,0.015),   # EC 1, changed to match Gard; 99th met, sed = 0.3
        Au=(0,0.06),    # EC 0.1, Gard has 0.06
        B=(0,0.1),      # 99th met = 0.88... 10x higher than everything else?
        Ba=(0,0.75),    # EC 1.0; 99th sed = 0.5
        Be=(0,0.007),   # EC 0.01; changed to match Gard
        Bi=(0,0.0007),  # EC 0.2; changed to match Gard. 99th sed = 0.041
        C=(0,10), 
        Cd=(0,0.001),   # EC 0.1, changed to match Gard
        Ce=(0,0.16),    # Gard has 0.1, but 99th sed = 0.094
        Cl=(0,8),       # EC 10; Gard has 0.7%; 99th ign = 2.5%, 99th sed = 6.4%
        Co=(0,1),       # EC 1; Gard has 0.025
        Cr2O3=(0,10),   # Pure magnesiochromite (MgCr₂O₄) is 79%, 99th ign = 2.1%
        Cs=(0,0.01),    # EC 0.009, changed to match Gard; 99th met = 0.011
        Cu=(0,0.5),     # EC 1; 99th sed = 0.48
        Dy=(0,0.007),   # EC 0.009; Gard has 0.007
        Er=(0,0.003),   # EC 0.0035; Gard has 0.0027; 99th sed = 0.0045. From histogram, 0.003 is reasonable
        Eu=(0,0.002),   # EC 0.0025; matched to Gard
        F=(0,1),        # EC ??; Gard has 0.5%; seds have more F (99th sed = 3.2%, 
                        # 99th ign = 0.77%), altho from the histogram, it's unclear how 
                        # much that's skewed by outliers. FWIW,  fluorite (CaF₂) is 25% F.
        Ga=(0,0.006),   # EC 0.01; 99th met = 0.0056
        Gd=(0,0.008),   # EC 0.012; Gard has 0.005; 99th sed = 0.0079
        Hf=(0,0.008),
        Hg=(0,0.0007),  # EC 0.01, changed to match Gard. Not really enough data
        Ho=(0,0.0011),  # EC 1e-6; Gard has 0.0011; 99th sed = 0.0019, 10x greater than others
        I=(0,0.001),    # Not enough data
        In=(0,0.002),   # EC 0.01; 99th sed = 0.0017
        Ir=(0,0.01),    # EC 0.001; Gard has 0.01; Not really enough data
        La=(0,0.1),     # EC 0.08, changed to match Gard
        Li=(0,0.1),     # EC 0.3; 99th sed = 0.083
        Lu=(0,0.0006),  # EC 0.0007; changed to match Gard 
        MnO=(0,4),      # EC 3; 99th sed = 3.3
        Mo=(0,0.2),     # EC 0.18, changed to match Gard; 99th met = 0.064
        Nb=(0,0.08),
        Nd=(0,0.04),    # EC 1e-5; Gard has 0.04; 99th sed = 0.031
        NiO=(0,10),     # EC 3, 99th met = 0.61, sed spike at ~0.65. Changed to match Gard
        Os=(0,0.01),    # Gard has 0.0001
        P2O5=(0,30),    # EC 40.7; pure fluorapatite (Ca₅(PO₄)₃F) is 28% P₂O₅; 99th sed = 18%. Unclear where EC max came from
        Pb=(0,0.04),    # EC 5, changed to match Gard. 99th met = 0.24
        Pd=(0,0.03),    # EC 0.01; Gard has 0.03
        Pt=(0,0.1),     # EC 0.01; Gard has 0.1
        Pr=(0,0.02),    # EC 1e-6; Gard has 0.01, 99th sed = 0.012
        Re=(0,1e-6),    # EC 0.01, Gard has 3.0e-7.... we divide anything smaller than 1e-6
        Rb=(0,0.12),
        Sb=(0,0.001),   # EC 0.01, changed to match Gard;
        Sc=(0,0.02),    # EC 0.01; changed to match Gard
        Se=(0,0.08),    # EC 0.1, Gard has 0.01; 99th sed = 0.05
        S=(0,28),       # EC 10; Gard has ~ 1.5 wt.%. Brenhin has 10%. 99th sed = 11%. Pyrite (FeS₂) = 27%
        Sm=(0,0.01),    # EC 2e-6; Gard has 0.01; 99th sed = 0.0085
        Sn=(0,0.005),   # EC 0.1, changed to match Gard; 99th sed / met = 0.02
        Sr=(0,0.7), Ta=(0,0.004),
        Tb=(0,0.0014),  # EC 0.0012; 99th sed = 0.0013
        Te=(0,0.01),
        Th=(0,0.02),    # EC 0.3, changed to match Gard; 99th sed = 0.013
        Tl=(0,0.0015),  # EC 0.1, changed to match Gard; 99th met = 0.24
        Tm=(0,0.0008),
        U=(0,0.02),     # EC 0.3, Gard has 0.005; 99th sed = 0.24; probably unit contam but set from histogram
        V=(0,0.15),
        W=(0,0.001),    # EC 0.1, Gard has 0.0005; we unit convert anything above 0.001
        Y=(0,0.045), Yb=(0,0.003), 
        Zn=(0,1),       # Gard has 0.07; 99th sed = 0.43
        Zr=(0,0.3),
    )

    # Restrict from bounds
    for k in keys(bounds)
        if haskey(bulk, k)
            # Replace 0 with NaN: a 0 is likely below detection and not *actually* 0
            bulk[k][.!(bounds[k][1] .< bulk[k] .<= bounds[k][2])] .= NaN
        elseif warn 
            @warn "Tuple bulk does not contain key $k."
        end
    end

    # Some bounds are based on element-oxide ratios or should have special restrictions
    # based on rock type
    include_minor!(cats)

    # Scandium, strontium, molybdenum, cesium
    bulk.Sc[bulk.Sc .* bulk.SiO2 .>= 5000] .= NaN
    bulk.Sr[bulk.Sr .* bulk.SiO2 .> 10^5.5] .= NaN
    bulk.Mo[bulk.Mo ./ bulk.SiO2.^4 .> 10^-4.2] .= NaN
    bulk.Cs[bulk.Cs ./ bulk.SiO2.^2 .> 0.03] .= NaN

    # Iron and titanium
    bulk.FeOT[(cats.ign .| cats.met) .& (bulk.FeOT .* bulk.SiO2 .> 1200)] .= NaN
    bulk.TiO2[(cats.ign .| cats.met) .& (bulk.TiO2 .* bulk.SiO2 .> 300)] .= NaN
    bulk.FeOT[cats.chert .& .!(0 .< bulk.FeOT .< 90)] .= NaN  # 99th for cherts (BIFs) is 89%
    bulk.FeOT[cats.sed .& .!(0 .< bulk.FeOT .< 30)] .= NaN
    bulk.FeOT[(cats.sed .| cats.met .| cats.ign) .& .!(0 .< bulk.FeOT .< 30)] .= NaN    # Other shit

    # Warn for non-screened keys 
    if warn 
        unscreened = "" 
        bk = keys(bounds)
        for k in keys(bulk)
            (!(k in bk) & (k != :FeOT)) && (unscreened *= "$k; ")
        end
        @warn "These keys were not screened for outliers: $unscreened."
    end

    return bulk
end
export screen_outliers!

## --- End of File