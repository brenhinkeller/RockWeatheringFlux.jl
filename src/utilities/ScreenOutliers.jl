"""
```julia
screen_outliers!(bulk)
```

Remove any geochemical outliers in `bulk` and replace with `NaN`s.

"""
function screen_outliers!(bulk::NamedTuple)
    # Oxide outliers
    bulk.SiO2[(0.001 .> bulk.SiO2) .| (bulk.SiO2 .> 100)] .= NaN
    bulk.TiO2[(0.01 .> bulk.TiO2) .| (bulk.TiO2 .> 10)] .= NaN      # bulk.TiO2.*bulk.SiO2.>300
    bulk.Al2O3[(0.1 .>= bulk.Al2O3) .| (bulk.Al2O3 .>= 30)] .= NaN
    bulk.MgO[(0.003 .>= bulk.MgO) .| (bulk.MgO .>= 50)] .= NaN      # Pure forsterite is 35#
    bulk.CaO[(0.003 .>= bulk.CaO) .| (bulk.CaO .>= 56)] .= NaN      # Max is pure CaCO3
    bulk.Na2O[(0.003 .>= bulk.Na2O) .| (bulk.Na2O .>= 30)] .= NaN
    bulk.K2O[(0.003 .>= bulk.K2O) .| (bulk.K2O .>= 30)] .= NaN
    bulk.P2O5[(0.001 .>= bulk.P2O5) .| (bulk.P2O5 .> 40.7)] .= NaN   # Max is pure CaPO4
    bulk.MnO[(0.001 .>= bulk.MnO) .| (bulk.MnO .> 3)] .= NaN
    bulk.Cr2O3[(0.0001 .>= bulk.Cr2O3) .| (bulk.Cr2O3 .>= 10)] .= NaN
    bulk.NiO[(0.0003 .>= bulk.NiO) .| (bulk.NiO .>= 3)] .= NaN
    bulk.FeOT[(bulk.FeOT .>= 0.01) .| (bulk.FeOT .> 100) .| (bulk.FeOT .* bulk.SiO2 .> 1200)] .= NaN
    
    # Element outliers
    bulk.Ag[(1e-2/10000 .>= bulk.Ag) .| (bulk.Ag .>= 3e3/10000)] .= NaN
    bulk.As[(0.1/10000 .>= bulk.As) .| (bulk.As .>= 1e4/10000)] .= NaN
    bulk.Au[(1e-4/10000 .>= bulk.Au) .| (bulk.Au .>= 1e3/10000)] .= NaN
    bulk.B[(0.1/10000 .>= bulk.B) .| (bulk.B .>= 1e3/10000)] .= NaN
    bulk.Ba[(1/10000 .>= bulk.Ba) .| (bulk.Ba .> 1e4/10000)] .= NaN
    bulk.Be[(0.03/10000 .>= bulk.Be) .| (bulk.Be .>= 100/10000)] .= NaN
    bulk.Bi[(1e-3/10000 .>= bulk.Bi) .| (bulk.Bi .>= 2e3/10000)] .= NaN
    bulk.C[(1e-3/10000 .>= bulk.C) .| (bulk.C .>= 1e5/10000)] .= NaN
    bulk.Cd[(1e-2/10000 .>= bulk.Cd) .| (bulk.Cd .>= 1e3/10000)] .= NaN
    
    bulk.Cl[(5/10000 .>= bulk.Cl) .| (bulk.Cl .>= 1e5/10000)] .= NaN
    bulk.Co[(0.1/10000 .>= bulk.Co) .| (bulk.Co .>= 1e4/10000)] .= NaN
    bulk.Cs[(0.001/10000 .> bulk.Cs) .| (bulk.Cs .>= 90/10000) .| (bulk.Cs ./ bulk.SiO2.^2 .> 2)] .= NaN
    bulk.Cu[(0.2/10000 .>= bulk.Cu) .| (bulk.Cu .>= 1e4/10000)] .= NaN
    bulk.Dy[(0.01/10000 .>= bulk.Dy) .| (bulk.Dy .>= 90/10000)] .= NaN
    bulk.Er[(0.02/10000 .>= bulk.Er) .| (bulk.Er .>= 35/10000)] .= NaN
    bulk.Eu[(0.01/10000 .>= bulk.Eu) .| (bulk.Eu .>= 25/10000)] .= NaN
    
    bulk.Ga[(1/10000 .>= bulk.Ga) .| (bulk.Ga .>= 1e2/10000)] .= NaN
    bulk.Gd[(0.01/10000 .>= bulk.Gd) .| (bulk.Gd .>= 120/10000)] .= NaN
    bulk.Hf[(2e-2/10000 .> bulk.Hf) .| (bulk.Hf .>= 80/10000)] .= NaN
    bulk.Hg[(1e-3/10000 .>= bulk.Hg) .| (bulk.Hg .>= 1e2/10000)] .= NaN
    
    bulk.I[(1e-4/10000 .>= bulk.I) .| (bulk.I .>= 1e2/10000)] .= NaN
    bulk.In[(1e-2/10000 .>= bulk.In) .| (bulk.In .>= 1e2/10000)] .= NaN
    bulk.Ir[(1e-6/10000 .> bulk.Ir) .| (bulk.Ir .>= 10/10000)] .= NaN
    bulk.La[(0.05/10000 .>= bulk.La) .| (bulk.La .>= 800/10000)] .= NaN
    bulk.Li[(0.3/10000 .>= bulk.Li) .| (bulk.Li .>= 3e3/10000)] .= NaN
    bulk.Lu[(0.005/10000 .>= bulk.Lu) .| (bulk.Lu .>= 07/10000)] .= NaN
    bulk.Mo[(0.05/10000 .>= bulk.Mo) .| (bulk.Mo .>= 1800/10000)] .= NaN
    bulk.Nb[(0.1/10000 .>= bulk.Nb) .| (bulk.Nb .> 800/10000)] .= NaN
    
    bulk.Os[(1e-8/10000 .> bulk.Os) .| (bulk.Os .>= 1e2/10000)] .= NaN
    bulk.Pb[(5e-2/10000 .>= bulk.Pb) .| (bulk.Pb .>= 5e4/10000)] .= NaN
    bulk.Pd[(1e-5/10000 .>= bulk.Pd) .| (bulk.Pd .>= 1e2/10000)] .= NaN
    bulk.Pt[(1e-5/10000 .>= bulk.Pt) .| (bulk.Pt .>= 1e2/10000)] .= NaN
    
    bulk.Re[(1e-8/10000 .> bulk.Re) .| (bulk.Re .>= 1e2/10000)] .= NaN
    bulk.Rb[(0.03/10000 .> bulk.Rb) .| (bulk.Rb .>= 1200/10000)] .= NaN
    bulk.Sb[(1e-2/10000 .>= bulk.Sb) .| (bulk.Sb .>= 1e2/10000)] .= NaN
    bulk.Sc[(0.2/10000 .>= bulk.Sc) .| (bulk.Sc .>= 100/10000)] .= NaN
    bulk.Se[(1e-2/10000 .>= bulk.Se) .| (bulk.Se .>= 1e3/10000)] .= NaN
    bulk.S[(1/10000 .>= bulk.S) .| (bulk.S .>= 1e5/10000)] .= NaN
    
    bulk.Sn[(0.1/10000 .>= bulk.Sn) .| (bulk.Sn .>= 1e3/10000)] .= NaN
    bulk.Sr[(1/10000 .>= bulk.Sr) .| (bulk.Sr .>= 7000/10000)] .= NaN
    bulk.Ta[(1e-3/10000 .>= bulk.Ta) .| (bulk.Ta .>= 40/10000)] .= NaN
    bulk.Tb[(0.01/10000 .>= bulk.Tb) .| (bulk.Tb .>= 12/10000)] .= NaN
    bulk.Te[(1e-4/10000 .>= bulk.Te) .| (bulk.Te .>= 1e2/10000)] .= NaN
    bulk.Th[(1e-3/10000 .>= bulk.Th) .| (bulk.Th .>= 3e3/10000)] .= NaN
    bulk.Tl[(1e-3/10000 .>= bulk.Tl) .| (bulk.Tl .>= 1e3/10000)] .= NaN
    bulk.Tm[(0.01/10000 .>= bulk.Tm) .| (bulk.Tm .>= 08/10000)] .= NaN
    bulk.U[(1e-4/10000 .>= bulk.U) .| (bulk.U .>= 3e3/10000)] .= NaN
    bulk.V[(1/10000 .>= bulk.V) .| (bulk.V .>= 1500/10000)] .= NaN
    bulk.W[(1e-2/10000 .>= bulk.W) .| (bulk.W .>= 1e3/10000)] .= NaN
    bulk.Y[(1/10000 .>= bulk.Y) .| (bulk.Y .> 450/10000)] .= NaN
    bulk.Yb[(0.02/10000 .>= bulk.Yb) .| (bulk.Yb .>= 30/10000)] .= NaN
    bulk.Zn[(1/10000 .>= bulk.Zn) .| (bulk.Zn .>= 1e4/10000)] .= NaN
    bulk.Zr[(1/10000 .>= bulk.Zr) .| (bulk.Zr .>= 3000/10000)] .= NaN

    # Why excluding oddly specific numbers?
    # bulk.Ce[bulk.Ce .>= 1600 .| bulk.Ce .<= 0.1] .= NaN; bulk.Ce[bulk.Ce==500] .= NaN
    # bulk.Ho[bulk.Ho .>= 15 .| bulk.Ho .<= 0.01] .= NaN; bulk.Ho[bulk.Ho==10] .= NaN
    # bulk.Nd[bulk.Nd .>= 400 .| bulk.Nd .<= 0.1] .= NaN ;bulk.Nd[bulk.Nd==300 .| bulk.Nd==200 .| bulk.Nd==100 .| bulk.Nd==150] .= NaN
    # bulk.Pr[bulk.Pr .>= 130 .| bulk.Pr .<= 0.01] .= NaN; bulk.Pr[bulk.Pr==68] .= NaN
    # bulk.Sm[bulk.Sm .>= 120 .| bulk.Sm .<= 0.02] .= NaN; bulk.Sm[bulk.Sm==50] .= NaN

    # All well and good, but 
    # bulk.F[bulk.F .>= 1e5 .| bulk.F .< 5] .= NaN # bimodal by four orders of magnitude [ppm vs wt.# splits between .<1 and .>20]
    

    return bulk
end

# bulk.Na[bulk.Na .> 1e5 .| bulk.Na .< 100] .= NaN
# bulk.Mg[bulk.Mg .> 2e5 .| bulk.Mg .< 15] .= NaN
# bulk.Al[bulk.Al .>= 2e5 .| bulk.Al .<= 1e2] .= NaN
# bulk.P[bulk.P .>= 2.4e5 .| bulk.P .<= 10] .= NaN # Max is pure CaPO4
# bulk.K[bulk.K .>= 1.5e5 .| bulk.K .<= 50] .= NaN
# bulk.Ca[bulk.Ca .>= 4e5 .| bulk.Ca .<= 50] .= NaN # Max is pure CaCO3 (and not CaO??)
# bulk.Ti[bulk.Ti .>= 3e4 .| bulk.Ti .<= 10] .= NaN
# bulk.Cr[bulk.Cr .>= 1e4 .| bulk.Cr .<= 0.3] .= NaN
# bulk.Mn[bulk.Mn .>= 3e5 .| bulk.Mn .<= 5] .= NaN #Max is pure MnO
# bulk.Fe[bulk.Fe .>= 7e5 .| bulk.Fe .<= 1e2] .= NaN #Max is pure Fe3O4
# bulk.Ni[bulk.Ni .>= 1e4 .| bulk.Ni .<= 0.2] .= NaN





















