## --- File names
    # 500 Macrostrat samples
        # macrostrat_io = "output/toy_responses.h5"
        # matchedbulk_io = "output/toy_bulkidx.tsv"

        # ucc_out = "results/toy_exposedcrust.tsv"
        # eroded_out = "output/toy_erodedmaterial.h5"
        # erodedabs_out = "results/toy_erodedmaterial_abs.tsv"
        # erodedrel_out = "results/toy_erodedmaterial_rel.tsv"

    # 50_000 Macrostrat samples (note: .tsv filetype is deprecated)
        # macrostrat_io = "output/pregenerated_responses.tsv"
        # matchedbulk_io = "output/bulkidx.tsv"

        # ucc_out = "results/exposedcrust.tsv"
        # eroded_out = "output/erodedmaterial.h5"
        # erodedabs_out = "results/erodedmaterial_abs.tsv"
        # erodedrel_out = "results/erodedmaterial_rel.tsv"

    # 250_000 Macrostrat samples
        macrostrat_io = "output/250K_responses.h5"
        matchedbulk_io = "output/250K_bulkidx.tsv"

        ucc_out = "results/250K_exposedcrust.tsv"
        eroded_out = "output/250K_erodedmaterial.h5"
        erodedabs_out = "results/250K_erodedmaterial_abs.tsv"
        erodedrel_out = "results/250K_erodedmaterial_rel.tsv"

    # 1_000_000 Macrostrat samples
        # macrostrat_io = "output/1M_responses.h5"
        # matchedbulk_io = "output/1M_bulkidx.tsv"

        # ucc_out = "results/1M_exposedcrust.tsv"
        # eroded_out = "output/1M_erodedmaterial.h5"
        # erodedabs_out = "results/1M_erodedmaterial_abs.tsv"
        # erodedrel_out = "results/1M_erodedmaterial_rel.tsv"


## --- Major and minor elements
    """
    ```julia
    get_elements()
    ```

    Define major and minor elements, returned as `Vector{Symbol}`s

    Major elements:
    * SiO2, Al2O3, FeOT, TiO2, MgO, CaO, Na2O, K2O

        
    Minor elements:
    * Ag, As, Au, B, Ba, Be, Bi, C, CaCO3, Cd, Ce, Cl, Co, Cr2O3, Cs, Cu, Dy, Er, Eu, 
        F, Ga, Gd, Hf, Hg, Ho, I, In, Ir, La, Li, Lu, MnO, Mo, Nb, Nd, NiO, Os, P2O5, Pb, 
        Pd, Pt, Pr, Re, Rb, Sb, Sc, Se, S, Sm, Sn, Sr, Ta, Tb, Te, Th, Tl, Tm, U, V, W, Y, 
        Yb, Zn, Zr

    Major elements are in part defined based on Faye and Ødegård 1975 
    (https://www.ngu.no/filearchive/NGUPublikasjoner/NGUnr_322_Bulletin_35_Faye_35_53.pdf).

    See also: `major_elements`

    """
    function get_elements()
        majors = [:SiO2,:Al2O3,:FeOT,:TiO2,:MgO,:CaO,:Na2O,:K2O,]
        minors = [:Ag,:As,:Au,:B,:Ba,:Be,:Bi,:C,:CaCO3,:Cd,:Ce,:Cl,:Co,:Cr2O3,:Cs,:Cu,
            :Dy,:Er,:Eu,:F,:Ga,:Gd,:Hf,:Hg,:Ho,:I,:In,:Ir,:La,:Li,:Lu,:MnO,:Mo,:Nb,:Nd,
            :NiO,:Os,:P2O5,:Pb,:Pd,:Pt,:Pr,:Re,:Rb,:Sb,:Sc,:Se,:S,:Sm,:Sn,:Sr,:Ta,:Tb,
            :Te,:Th,:Tl,:Tm,:U,:V,:W,:Y,:Yb,:Zn,:Zr
        ]

        return majors, minors
    end


## --- Igneous rock classifications by silica content
    """
    ```
    get_ignsilica()
    ```

    Define felsic, intermediate, and mafic rocks by wt.% silica, using values from Keller 
    and Schoene, 2012 (DOI: 10.1038/nature11024):
      * Felsic: 62-74%
      * Intermediate: 51-62%
      * Mafic: 43-51%
    """
    get_ignsilica() = return (fel = (62, 74), int = (51, 62), maf = (43, 51))


## --- Rock type classifications based on EarthChem system in RockNameInference.m
    """
    ```julia
    get_rock_class(major::Bool, npoints::Int64)
    ```
    
    Define sedimentary, igneous, and metamorphic rock types and subtypes, and initialize a
    NamedTuple of `npoints`-element BitVectors for each defined rock type.
    
    If `major` is `true`, only major types are defined. If `major` is `false`, all subtypes
    are defined.

    See also: `get_minor_types`

    # Example
    ```julia
    typelist, cats = get_rock_class(true, 10)
    ```
    """
    function get_rock_class(major::Bool, npoints::Int64)
        # Sedimentary
        siliciclast = ("siliciclast", "conglomerat", "sand", "psamm", "arenit", "arkos", 
            "silt")
        shale = ("mud", "clay","shale", "wacke", "argillite", "argillaceous", "flysch", 
            "pelit", "turbidite")
        carb = ("carbonate", "limestone", "dolo", "marl", "chalk", "travertine", "tavertine", 
            "teravertine", "tufa")
        chert = ("chert", "banded iron")
        evaporite = ("evaporite", "gypsum", "salt flat", "caliche")
        coal = ("coal", "anthracite")

        sed = (("sediment", "fluv", "clast", "gravel", "pebble", "boulder", "diamict",
            "tillite", "stream", "beach", "terrace",  "marine deposits",  "paleosol")...,
            siliciclast..., shale..., carb..., chert..., evaporite..., coal...
        )

        # Igneous
        volc = ("volcan", "lava", "lahar", "ignimbrite", "ashfall", "tuff", "diatreme",
            "pipe", "basalt", "andesit", "dacit", "rhyolit", "pillow", "carbonatite", 
            "tephra", "obsidian", "ash", "scoria", "pumice", "cinder", "latite", 
            "basanite", "phonolite", "fonolito", "trachyte", "palagonite", "mugearite", 
            "kimberlite", "ultramafitite", "komatiite",)
        plut = (
            # "True" plutonic rocks
            ("pluton", "batholith", "granit", "tonalit", "gabbro", "norite", 
            "diorit", "monzonit", "syenit", "peridot", "dunit", "harzburg", "anorthosite", 
            "mangerite", "charnockite", "pegmatite", "aplite", "trond", "essexite", 
            "pyroxenite", "adamellite", "porphyry", "megacryst", "rapakivi", "bronzitite", 
            "alaskite", "troctolite")..., 
            # Hypabyssal rocks
            ("intrus", "hypabyssal", "sill", "dike", "stock", "laccolith", "lopolith", 
            "dolerit", "diabase", "porphyry", "microgranite")...
        )

        ign = (("igneous", "silicic ", "mafic", "felsic", "basite",)..., volc..., plut...)

        # Metamorphic
        metased = ("para", "metased", "meta-sed", "quartzite", "marble", "slate", "phyllite",)
        metaign = ("ortho", "metaign", "meta-ign", "serpentin", "amphibolit", "greenstone", 
            "eclogite", "metabasite",)
        lowgrade = ("slate", "phyllite", "serpentin", "greenstone", "greenschist", "zeolite", 
            "gossan", "alter", "hydrothermal", "palagonite",)
        highgrade = ("crystalline", "basement", "marble", "skarn", "schist", "blueschist", "gneiss", 
            "amphibolit", "eclogite", "granulit", "hornfels", "granofels", "sanidinite", "migma", 
            "enderbite", "anorthosite", "charnockite", "pyroxenite", "peridot", "dunit", "harzburg", 
            "high grade metamorphic")
        cataclastic = ("mylonit", "cataclasite", "melange", "gouge", "tecton",)

        met = (cataclastic..., ("meta", "calc silicate",)...,  metased..., metaign..., lowgrade..., 
            highgrade..., cataclastic...
        )

        # Cover
        cover = ("cover", "unconsolidated", "quaternary", "lluv", "soil", "regolith", 
            "laterite", "surficial deposits", "talus", "scree", "mass-wasting", "slide", 
            "peat", "swamp", "marsh", "water", "ice", "glaci", "till", "loess", "gravel", 
            "debris"
        )

        # Initialize type lists and BitVectors
        if major
            typelist = (sed=sed, ign=ign, met=met, cover=cover)
        else
            typelist = (siliciclast=siliciclast, shale=shale, carb=carb, chert=chert, 
                evaporite=evaporite, coal=coal, sed=sed, volc=volc, plut=plut, ign=ign, 
                metased=metased, metaign=metaign, met=met, cover=cover)
        end

        return typelist, NamedTuple{keys(typelist)}([falses(npoints) for _ in 1:length(typelist)]) 
    end

## --- Define minor types
    """
    ```julia
    get_minor_types()
    ```

    Return types nested under the sed, ign, and met "major" types.

    ### Minor Types:
      * Sed: siliciclast, shale, carb, chert, evaporite, coal
      * Ign: volc, plut
      * Met: metased, metaign
    
    # Example
    ```julia
    minorsed, minorign, minormet = get_minor_types()
    ```
    """
    get_minor_types() = (:siliciclast, :shale, :carb, :chert, :evaporite, :coal), 
        (:volc, :plut), (:metased, :metaign)


## --- Macrostrat type exclusions to avoid multi-matching
    """
    ```julia
    un_multimatch!(cats, major::Bool)
    ```

    Exclude Macrostrat matches from each other so each sample is only classified as one
    rock type.

    If `major` is `true`: 
    * Cover is excluded from all rock types.
    * Rocks classified as both metamorphic and sedimentary / igneous (i.e., metasedimentary
        and metaigneous rocks) are respectively re-classified as sedimentary and igneous rocks.
    * Rocks classified as sedimentary and igneous are re-classified as igneous rocks.

    If `major` is false:
    * Cover is excluded from all rock types.
    * Rocks classified as both metamorphic and sedimentary / igneous (i.e., metasedimentary
        and metaigneous rocks) are **excluded** from sedimentary and igneous rocks, and 
        re-classified as metasedimentary and metaigneous rocks.
    * Rocks classified as more than one subtype of sedimentary rocks are excluded from
        each other, in arbitrary order.
    * Rocks classified as both volcanic and plutonic, or both metasedimentary and 
        metaigneous, are respectively re-classified as undifferentiated igneous and metamorphic
        rocks.
    * Rocks classified as sedimentary and igneous are re-classified as igneous rocks.
    """
    un_multimatch!(cats, major::Bool) = _un_multimatch!(cats, static(major))
    function _un_multimatch!(cats, major::True)
        # Exclude cover
        cats.sed .&= .! cats.cover
        cats.ign .&= .! cats.cover
        cats.met .&= .! cats.cover

        # Classify metased as sed and metaign and ign
        cats.met .&= .! cats.sed
        cats.ign .&= .! cats.ign

        # Sed / ign rocks are classified as ign
        cats.sed .&= .! cats.ign

        return cats
    end
    
    function _un_multimatch!(cats, major::False)
        # Define types
        minorsed, minorign, minormet = get_minor_types()
        minortypes = (minorsed..., minorign..., minormet...)

        # Exclude cover from all major and minor rock types
        cats.sed .&= .! cats.cover
        cats.ign .&= .! cats.cover
        cats.met .&= .! cats.cover
        for i in minortypes
            cats[i] .&= .! cats.cover
        end

        # Exclude metamorphic rocks from sed and igns. Class as metased and metaign
        cats.metased .|= (cats.sed .& cats.met)
        cats.metaign .|= (cats.ign .& cats.met)

        cats.sed .&= .! cats.met
        for i in minorsed
            cats[i] .&= .! cats.met
        end

        cats.ign .&= .! cats.met
        for i in minorign
            cats[i] .&= .! cats.met
        end

        # Exclude sed subtypes from other sed subtypes
        for i in 1:(length(minorsed)-1)
            for j in minorsed[i+1:end]
                cats[minorsed[i]] .&= .! cats[j]
            end
        end

        # Both volcanic and plutonic is undifferentiated igneous
        cats.volc .&= .!(cats.volc .& cats.plut)
        cats.plut .&= .!(cats.volc .& cats.plut)

        # Both metaigneous and metasedimentary is undifferentiated metamorphic
        cats.metased .&= .!(cats.metased .& cats.metaign)
        cats.metaign .&= .!(cats.metased .& cats.metaign)

        # Sed / ign rocks are classified as ign
        cats.sed .&= .! cats.ign            
        for i in minorsed
            cats[i] .&= .! cats.ign 
        end

        return cats
    end

## --- End of file