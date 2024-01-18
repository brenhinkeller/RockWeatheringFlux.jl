## --- Create a test dataset 
    majors, minors = get_elements()
    allelements = [majors; minors]

    bulk = NamedTuple{Tuple(allelements)}(rand(-1:0.000001:120, 10) for _ in allelements)
    bulk_initial = deepcopy(bulk)


## --- Screen outliers and check that the function mutated the structure
    screen_outliers!(bulk)
    @test bulk_initial != bulk

## --- Test element oxides
    # I guess?
    tmin = 0
    tmax = 100

    @test all(isnan.(bulk.SiO2[.!(tmin .< bulk_initial.SiO2 .< tmax)]))
    @test all(isnan.(bulk.Al2O3[.!(tmin .< bulk_initial.Al2O3 .< tmax)]))
    @test all(isnan.(bulk.FeOT[.!(tmin .< bulk_initial.FeOT .< tmax)]))
    @test all(isnan.(bulk.TiO2[.!(tmin .< bulk_initial.TiO2 .< tmax)]))
    @test all(isnan.(bulk.MgO[.!(tmin .< bulk_initial.MgO .< tmax)]))
    @test all(isnan.(bulk.CaO[.!(tmin .< bulk_initial.CaO .< tmax)]))
    @test all(isnan.(bulk.Na2O[.!(tmin .< bulk_initial.Na2O .< tmax)]))
    @test all(isnan.(bulk.K2O[.!(tmin .< bulk_initial.K2O .< tmax)]))

    # Trace element oxides
    @test all(isnan.(bulk.P2O5[.!(tmin .< bulk_initial.P2O5 .< tmax)]))
    @test all(isnan.(bulk.MnO[.!(tmin .< bulk_initial.MnO .< tmax)]))
    @test all(isnan.(bulk.Cr2O3[.!(tmin .< bulk_initial.Cr2O3 .< tmax)]))
    @test all(isnan.(bulk.NiO[.!(tmin .< bulk_initial.NiO .< tmax)]))

##--- Test trace elements
    tmax = 10       # Maximum allowable trace concentration is 10 wt.%
    tmin = 1e-10    # Minimum allowable trace concentration is 1e-10 wt.%

    # for m in minors
    #     println("@test all(isnan.(bulk.$m[.!(tmin .< bulk_initial.$m .< tmax)]))")
    # end

    @test all(isnan.(bulk.Ag[.!(tmin .< bulk_initial.Ag .< tmax)]))
    @test all(isnan.(bulk.As[.!(tmin .< bulk_initial.As .< tmax)]))
    @test all(isnan.(bulk.Au[.!(tmin .< bulk_initial.Au .< tmax)]))
    @test all(isnan.(bulk.B[.!(tmin .< bulk_initial.B .< tmax)]))
    @test all(isnan.(bulk.Ba[.!(tmin .< bulk_initial.Ba .< tmax)]))
    @test all(isnan.(bulk.Be[.!(tmin .< bulk_initial.Be .< tmax)]))
    @test all(isnan.(bulk.Bi[.!(tmin .< bulk_initial.Bi .< tmax)]))
    @test all(isnan.(bulk.C[.!(tmin .< bulk_initial.C .< tmax)]))
    @test all(isnan.(bulk.Cd[.!(tmin .< bulk_initial.Cd .< tmax)]))
    @test all(isnan.(bulk.Ce[.!(tmin .< bulk_initial.Ce .< tmax)])) broken=true
    @test all(isnan.(bulk.Cl[.!(tmin .< bulk_initial.Cl .< tmax)]))
    @test all(isnan.(bulk.Co[.!(tmin .< bulk_initial.Co .< tmax)]))
    @test all(isnan.(bulk.Cs[.!(tmin .< bulk_initial.Cs .< tmax)]))
    @test all(isnan.(bulk.Cu[.!(tmin .< bulk_initial.Cu .< tmax)]))
    @test all(isnan.(bulk.Dy[.!(tmin .< bulk_initial.Dy .< tmax)]))
    @test all(isnan.(bulk.Er[.!(tmin .< bulk_initial.Er .< tmax)]))
    @test all(isnan.(bulk.Eu[.!(tmin .< bulk_initial.Eu .< tmax)]))
    @test all(isnan.(bulk.F[.!(tmin .< bulk_initial.F .< tmax)])) broken=true
    @test all(isnan.(bulk.Ga[.!(tmin .< bulk_initial.Ga .< tmax)]))
    @test all(isnan.(bulk.Gd[.!(tmin .< bulk_initial.Gd .< tmax)]))
    @test all(isnan.(bulk.Hf[.!(tmin .< bulk_initial.Hf .< tmax)]))
    @test all(isnan.(bulk.Hg[.!(tmin .< bulk_initial.Hg .< tmax)]))
    @test all(isnan.(bulk.Ho[.!(tmin .< bulk_initial.Ho .< tmax)])) broken=true
    @test all(isnan.(bulk.I[.!(tmin .< bulk_initial.I .< tmax)]))
    @test all(isnan.(bulk.In[.!(tmin .< bulk_initial.In .< tmax)]))
    @test all(isnan.(bulk.Ir[.!(tmin .< bulk_initial.Ir .< tmax)]))
    @test all(isnan.(bulk.La[.!(tmin .< bulk_initial.La .< tmax)]))
    @test all(isnan.(bulk.Li[.!(tmin .< bulk_initial.Li .< tmax)]))
    @test all(isnan.(bulk.Lu[.!(tmin .< bulk_initial.Lu .< tmax)]))
    @test all(isnan.(bulk.Mo[.!(tmin .< bulk_initial.Mo .< tmax)]))
    @test all(isnan.(bulk.Nb[.!(tmin .< bulk_initial.Nb .< tmax)]))
    @test all(isnan.(bulk.Nd[.!(tmin .< bulk_initial.Nd .< tmax)])) broken=true
    @test all(isnan.(bulk.Os[.!(tmin .< bulk_initial.Os .< tmax)]))
    @test all(isnan.(bulk.Pb[.!(tmin .< bulk_initial.Pb .< tmax)]))
    @test all(isnan.(bulk.Pd[.!(tmin .< bulk_initial.Pd .< tmax)]))
    @test all(isnan.(bulk.Pt[.!(tmin .< bulk_initial.Pt .< tmax)]))
    @test all(isnan.(bulk.Pr[.!(tmin .< bulk_initial.Pr .< tmax)])) broken=true
    @test all(isnan.(bulk.Re[.!(tmin .< bulk_initial.Re .< tmax)]))
    @test all(isnan.(bulk.Rb[.!(tmin .< bulk_initial.Rb .< tmax)]))
    @test all(isnan.(bulk.Sb[.!(tmin .< bulk_initial.Sb .< tmax)]))
    @test all(isnan.(bulk.Sc[.!(tmin .< bulk_initial.Sc .< tmax)]))
    @test all(isnan.(bulk.Se[.!(tmin .< bulk_initial.Se .< tmax)]))
    @test all(isnan.(bulk.S[.!(tmin .< bulk_initial.S .< tmax)]))
    @test all(isnan.(bulk.Sm[.!(tmin .< bulk_initial.Sm .< tmax)])) broken=true
    @test all(isnan.(bulk.Sn[.!(tmin .< bulk_initial.Sn .< tmax)]))
    @test all(isnan.(bulk.Sr[.!(tmin .< bulk_initial.Sr .< tmax)]))
    @test all(isnan.(bulk.Ta[.!(tmin .< bulk_initial.Ta .< tmax)]))
    @test all(isnan.(bulk.Tb[.!(tmin .< bulk_initial.Tb .< tmax)]))
    @test all(isnan.(bulk.Te[.!(tmin .< bulk_initial.Te .< tmax)]))
    @test all(isnan.(bulk.Th[.!(tmin .< bulk_initial.Th .< tmax)]))
    @test all(isnan.(bulk.Tl[.!(tmin .< bulk_initial.Tl .< tmax)]))
    @test all(isnan.(bulk.Tm[.!(tmin .< bulk_initial.Tm .< tmax)]))
    @test all(isnan.(bulk.U[.!(tmin .< bulk_initial.U .< tmax)]))
    @test all(isnan.(bulk.V[.!(tmin .< bulk_initial.V .< tmax)]))
    @test all(isnan.(bulk.W[.!(tmin .< bulk_initial.W .< tmax)]))
    @test all(isnan.(bulk.Y[.!(tmin .< bulk_initial.Y .< tmax)]))
    @test all(isnan.(bulk.Yb[.!(tmin .< bulk_initial.Yb .< tmax)]))
    @test all(isnan.(bulk.Zn[.!(tmin .< bulk_initial.Zn .< tmax)]))
    @test all(isnan.(bulk.Zr[.!(tmin .< bulk_initial.Zr .< tmax)]))


## --- End of file 