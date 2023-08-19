# Definitions, file names, etc. Contains most hard-coded values.

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


## --- Color names for visualization

    clr_gradient = :jet1

    # Observed samples by rock type
    clr_ign = :crimson
    clr_volc = :tomato 
    clr_plut = :darkred
    clr_sed = :green

    # Resampled
    clr_rs = :grey


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


## --- Rock type classifications
    """
    ```julia
    get_rock_class([major::Bool], [inclusive::Bool])
    ```

    Define sedimentary, igneous, and metamorphic rock types and subtypes.

    ### Optional kwargs
        
        major
    
    Return only Tuples for sedimentary, metamorphic, and igneous types. Boolean; defaults 
    to `false`

        inclusive
    
    Sedimentary, igneous, and metamorphic lists include terms for all included subtypes. 
    Boolean; defaults to `false` unless `major` is `true`.
    """
    function get_rock_class(major::Bool=false, inclusive::Bool=false)
        # Sedimentary
        siliciclast = ("siliciclast", "conglo", "sand", "psamm", "arenit", "arkos", "silt")
        shale = ("lutite", "mud", "clay", "shale", "wacke", "argillite", "argillaceous", 
            "flysch", "pelit", "turbidite", "tasmanite", "breccia")
        carb = ("carbonate", "limestone", "dolo", "marl", "chalk", "coquina", "biogenic", 
            "travertine", "tavertine", "tufa", "calcarenite", "teravertine")
        chert = ("chert", "opal", "porcellanite", "diatomite", "novaculite", "iron", "taconite", 
            "banded iron")
        evaporite = ("evaporite", "anhydrite", "gypsum", "trona", "halite", "sylvite", 
            "salt flat", "caliche")
        phosphorite = ("phosphorite", "phosphate")
        coal = ("coal", "anthracite", "peat", "lignite", "bitumen")
        volcaniclast = ("tonstein", "peperite", "volcaniclastic")
        sed = ("sediment", "clast", "diamict","tillite", "stream  deposits", "beach deposits", 
            "terrace",  "marine deposits",  "paleosol", "spiculite", "glauconite")

        # Igneous
        volc = ("volcanic", "extrusive", "tuff", "basalt", "andesit", "dacit", "rhyolit", 
            "pillow", "glass", "tephra", "obsidian", "ash", "scoria", "pumice", "cinder", 
            "lahar", "lava", "lenticulite", "absarokite", "adakite", "alnoite", "alvikite", 
            "analcimite", "anamesite", "ankaramite", "ankaratrite", "augitite", "basanit", 
            "benmoreite", "bergalite", "boninite", "breunneritite", "buchonite", "campanite", 
            "camptonite", "ciminite", "damkjernite", "dellenite", "domite", "etindite", 
            "eutaxite", "felsite", "foidite", "fortunite", "gauteite", "grazinite", "hauynophyre", 
            "hawaiite", "icelandite", "ignimbrite", "kalsilit", "katungite", "kenyte", 
            "keratophyre", "kersantite", "kimberlite", "kivite", "komatiite", "lampro", "latite", 
            "leucitite", "liparite", "limburgite", "linosaite", "madupite", "mafite", "mafraite", 
            "mafurite", "marianite", "meimechite", "melafoidite", "melilit", "melnoite", 
            "minette", "monchiquite", "mondhaldeite", "mugearite", "murambite", "nephelinite", 
            "oceanite", "orangeite", "ordanchite", "orendite", "picrite", "phonolite", 
            "pantellerite", "polzenite", "palagonite", "sannaite", "santorinite", "shoshonite", 
            "tachylyte", "tahitite", "tephrite", "tholeiite", "trachyte", "tristanite", 
            "vicoite", "vitrophere", "vitrophyre", "vulsinite", "wyomingite", "volcan", 
            "ashfall", "diatreme", "pipe", "carbonatite", "ash fall", "basanite", "fonolito", 
            "ultramafitite")
        plut = ("plutonic", "intrusive", "granit", "tonalit", "gabbro", "diorit", "monzonit", 
            "syenit", "adamellite", "alaskite", "allivalite", "anorthosite", "apatitite", 
            "aplite", "biotitite", "borolanite", "bostonite", "bronzitite", "chromitite", 
            "comendite", "corganite", "corgaspinite", "cortlandite", "crinanite", "charnockite", 
            "diabase", "dolerit", "dunit", "durbachite", "enderbite", "essexite", "fergusite", 
            "foyaite", "glenmuirite", "glimmerite", "granophyre", "hypabyssal", "harzburg", 
            "hauynite", "hornblendite", "ijolite", "jacupirangite", "jotunite", "juvite", 
            "kaersutitite", "kamafugite", "kentallenite", "kullaite", "labradorite", "larvikite", 
            "lestiwarite", "lherzolite", "lujavrite", "luscladite", "marscoite", "melteigite", 
            "megacryst", "missourite", "norite", "nordmarkite", "olivinite", "orthoclasite", 
            "ottajanite", "pegmat", "peridot", "pyroxenite", "porphyry", "porphyrite", "puglianite",
            "riedenite", "sanukite", "sebastianite", "shonkinite", "sommaite", "sovite", 
            "tannbuschite", "teschenite", "theralite", "tinguaite", "trond", "vsbergite", 
            "topazite", "troctolite", "turjaite", "ugandite", "uncompahgrite", "urtite", 
            "vaugnerite", "vibetoite", "websterite", "wehrlite", "yamaskite", "pluton", 
            "batholith", "mangerite", "pegmatite", "rapakivi", "intrus", "sill", "dike", 
            "stock", "laccolith", "lopolith", "microgranite")
        ign = ("igneous", "silicic ", "mafic", "felsic", "basite", "phoscorite", "rauhaugite", 
            "beforsite")

        # Metamorphic
        metased = ("para", "metased", "quartzite", "marble", "slate", "leptite", "phyllite", 
            "porcellanite", "meta-sed", "hornfels")
        metaign = ("orthogneiss", "metaign", "serpentin", "amphibolit", "greenstone", "eclogite", 
            "basite", "greisen", "halleflinta", "leucophyre", "melaphyre", "propylite", "spilite", 
            "ultramafitite", "alkremite", "ortho", "meta-ign", "metabasite")
        lowgrade = ("slate", "phyllite", "serpentin", "greenstone", "greenschist", "zeolite", 
            "gossan", "alter", "hydrothermal", "palagonite",)
        highgrade = ("crystalline", "basement", "marble", "skarn", "schist", "blueschist", 
            "gneiss", "amphibolit", "eclogite", "granulit", "granofels", "sanidinite", 
            "migma", "enderbite", "anorthosite", "charnockite", "pyroxenite", "peridot", "dunit", 
            "harzburg", "high grade metamorphic")
        cataclastic = ("mylonit", "cataclasite", "melange", "gouge", "tecton",)
        met = (("meta", "garnet", "buchite", "epidot", "fenite", "albitite", "chloritite", 
            "phlogopitite", "calc silicate", "calcsilicate", "rodingite", "sericitite", 
            "tactite", "soapstone", "talc", "tourmalinite", "unakite", "vogesite")...,
            lowgrade..., highgrade..., cataclastic...,
        )

        # Cover
        cover = ("lluv", "fluv", "boulder", "gravel", "aleurite", "glaci", "till", "loess", 
            "regolith", "debris", "fill", "slide", "unconsolidated", "talus", "stream", "beach",
            "terrace", "placer", "paleosol", "mass-wasting", "pebble", "cover", "quaternary", 
            "soil", "laterite", "surficial deposits", "scree", "peat", "swamp", "marsh", 
            "water", "ice")

        # If inclusive or major classes only, include subtypes
        if inclusive || major
            sed = unique((sed..., siliciclast..., shale..., carb..., chert..., evaporite..., 
                coal..., phosphorite..., volcaniclast...,))
            ign = unique((ign..., volc..., plut...))
            met = unique((met..., metased..., metaign..., lowgrade..., highgrade..., 
                cataclastic...))

            # If major elements only, return major elements
            major && return (sed=sed, ign=ign, met=met, cover=cover)
        end

        # Otherwise, return Tuple where major types do not include subtypes
        return (siliciclast = siliciclast, shale = shale, carb = carb, chert = chert, 
            evaporite = evaporite, coal = coal, phosphorite = phosphorite, 
            volcaniclast = volcaniclast, sed = sed, 
            volc = volc, plut = plut, ign = ign, 
            metased = metased, metaign = metaign, met = met, 
            cover = cover
        )
    end

    """
    ```julia
    get_cats(major::Bool, npoints::Int64)
    ```
    
    Initialize a NamedTuple of `npoints`-element BitVectors for each defined rock type.
    
    If `major` is `true`, only major types are defined. If `major` is `false`, all subtypes
    are defined.

    See also: `get_minor_types`

    # Example
    ```julia
    typelist, cats = get_cats(true, 10)
    ```
    """
    function get_cats(major::Bool, npoints::Int64)
        typelist = get_rock_class(major)

        return typelist, NamedTuple{keys(typelist)}([falses(npoints) for _ in 1:length(typelist)]) 
    end

## --- Define minor types

    """
    ```julia
    get_minor_types()
    ```

    Return types nested under the sed, ign, and met "major" types.

    **Important: this function will break if major types are listed before minor types in
    the `get_cats` function, and if sedimentary types are not listed first!**

    ### Minor Types:
      * Sed: siliciclast, shale, carb, chert, evaporite, coal, phosphorite, volcaniclast
      * Ign: volc, plut
      * Met: metased, metaign
    
    # Example
    ```julia
    minorsed, minorign, minormet = get_minor_types()
    ```
    """
    function get_minor_types()
        types = get_cats(false, 1)[2]
        allkeys = collect(keys(types))

        sed = findfirst(==(:sed), allkeys)
        ign = findfirst(==(:ign), allkeys)
        met = findfirst(==(:met), allkeys)

        return Tuple(allkeys[1:sed-1]), Tuple(allkeys[sed+1:ign-1]), Tuple(allkeys[ign+1:met-1])
    end


## --- Rock type exclusions to avoid multi-matching

    """
    ```julia
    un_multimatch!(cats, major::Bool)
    ```

    Exclude rock type matches from each other so each sample is only classified as one
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