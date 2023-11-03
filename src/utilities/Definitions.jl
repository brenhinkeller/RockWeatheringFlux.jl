# Definitions, file names, etc. Contains most hard-coded values.

## --- File names

    # 500 Macrostrat samples
        # macrostrat_io = "output/toy_responses.h5"
        # matchedbulk_io = "output/toy_bulkidx.tsv"

        # ucc_out = "results/toy_exposedcrust.tsv"
        # eroded_out = "output/toy_erodedmaterial.h5"
        # erodedabs_out = "results/toy_erodedmaterial_abs.csv"
        # erodedrel_out = "results/toy_erodedmaterial_rel.csv"

    # 50_000 Macrostrat samples (note: .tsv filetype is deprecated)
        # macrostrat_io = "output/pregenerated_responses.tsv"
        # matchedbulk_io = "output/bulkidx.tsv"

        # ucc_out = "results/exposedcrust.tsv"
        # eroded_out = "output/erodedmaterial.h5"
        # erodedabs_out = "results/erodedmaterial_abs.csv"
        # erodedrel_out = "results/erodedmaterial_rel.csv"

    # 250_000 Macrostrat samples
        macrostrat_io = "output/250K_responses.h5"
        matchedbulk_io = "output/250K_bulkidx.tsv"

        ucc_out = "results/250K_exposedcrust.tsv"
        eroded_out = "output/250K_erodedmaterial.h5"
        erodedabs_out = "results/250K_erodedmaterial_abs.csv"
        erodedrel_out = "results/250K_erodedmaterial_rel.csv"

    # 1_000_000 Macrostrat samples
        # macrostrat_io = "output/1M_responses.h5"
        # matchedbulk_io = "output/1M_bulkidx.tsv"

        # ucc_out = "results/1M_exposedcrust.tsv"
        # eroded_out = "output/1M_erodedmaterial.h5"
        # erodedabs_out = "results/1M_erodedmaterial_abs.csv"
        # erodedrel_out = "results/1M_erodedmaterial_rel.csv"


## --- Color names for visualization

    c_gradient = :jet1

    # Observed samples by rock type
    colors = (
        ign = :crimson, volc = :tomato, plut = :deeppink,
        sed = :royalblue, siliciclast = :dodgerblue, shale = :cadetblue, carb = :aqua,
            chert = :deepskyblue, evaporite = :lightskyblue, coal = :midnightblue, 
            phosphorite = :turquoise, volcaniclast = :rebeccapurple,
        met = :peru, metased = :yellowgreen, metaign = :darkorange,
        cover = :grey
    )

    # Resampled
    c_rs = :grey


## --- Major and minor elements

    """
    ```julia
    get_elements()
    ```

    Define major and minor elements, returned as `Vector{Symbol}`s

    Major elements:
    * SiO₂, Al₂O₃, FeOT, TiO₂, MgO, CaO, Na₂O, K₂O, Volatiles

        
    Minor elements:
    * Ag, As, Au, B, Ba, Be, Bi, C, Cd, Ce, Cl, Co, Cr₂O₃, Cs, Cu, Dy, Er, Eu, 
        F, Ga, Gd, Hf, Hg, Ho, I, In, Ir, La, Li, Lu, MnO, Mo, Nb, Nd, NiO, Os, P₂O₅, Pb, 
        Pd, Pt, Pr, Re, Rb, Sb, Sc, Se, S, Sm, Sn, Sr, Ta, Tb, Te, Th, Tl, Tm, U, V, W, Y, 
        Yb, Zn, Zr

    Major elements are in part defined based on Faye and Ødegård 1975 
    (https://www.ngu.no/filearchive/NGUPublikasjoner/NGUnr_322_Bulletin_35_Faye_35_53.pdf).

    See also: `major_elements`

    """
    function get_elements()
        majors = [:SiO2,:Al2O3,:FeOT,:TiO2,:MgO,:CaO,:Na2O,:K2O,:Volatiles]
        minors = [:Ag,:As,:Au,:B,:Ba,:Be,:Bi,:C,:Cd,:Ce,:Cl,:Co,:Cr2O3,:Cs,:Cu,
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

    Define sedimentary, igneous, and metamorphic rock types and subtypes. Metamorphic rocks
    are grouped with their protoliths when possible.

    ### Possible subtypes
     *  Sedimentary: siliciclastic, shale, carbonate, evaporite, chert, phosphorite, coal,
        sedimentary (uncategorized)
     *  Igneous: 
         *  Volcanic: komatiite, basalt, andesite, dacite, rhyolite, alkaline, 
            volcaniclastic, volcanic (uncategorized)
         *  Plutonic: peridotite, pyroxenite, gabbro, diorite, trondhjemite, tonalite, 
            granodiorite, granite, alkaline, plutonic (uncategorized)
         *  Carbonatite
         *  Igneous (uncategorized)
      * Metamorphic: metamorphic (uncategorized)

    ### Optional kwargs
        
        major
    
    Return only Tuples for sedimentary, metamorphic, and igneous types. Boolean; defaults 
    to `false`

        inclusive
    
    Sedimentary, igneous, and metamorphic lists include terms for all included subtypes. 
    Boolean; defaults to `false` unless `major` is `true`.

    Note that volcanic and plutonic rocks will include all nested volcanic / plutonic 
    subtypes.
    """
    function get_rock_class(; major::Bool=false, inclusive::Bool=false)
        # Sedimentary
        siliciclast = ("siliciclast", "conglo", "sand", "psamm", "arenit", "arkos", "silt",
            "breccia", "quartzite")
        shale = ("lutite", "mud", "clay", "shale", "wacke", "argillite", "argillaceous", 
            "flysch", "pelit", "turbidite", "tasmanite", "slate", "phyllite",)
        carb = ("carbonate", "limestone", "dolo", "marl", "chalk", "coquina", "biogenic", 
            "travertine", "tavertine", "tufa", "calcarenite", "teravertine", "marble", 
            "calc silicate", "calcsilicate", "skarn", )
        evap = ("evaporite", "anhydrite", "gypsum", "trona", "halite", "sylvite", 
            "salt flat", "caliche", "exhalite")
        chert = ("chert", "opal", "porcellanite", "diatomite", "novaculite", "iron", 
            "taconite", "banded iron")
        phosphorite = ("phosphorite", "phosphate")
        coal = ("coal", "anthracite", "lignite", "bitumen")
        sed = ("sediment", "clast", "diamict", "tillite", "stream deposits", 
            "beach deposits", "terrace",  "marine deposits",  "paleosol", "spiculite", 
            "glauconite", "meta-sed", "metased", "paragneiss", "para")

        # Volcanic
        komatiite = ("komatiite", "meimechite", "ultramafitite")
        basalt = ("basalt", "pillow", "scoria", "picrite", "anamesite", "hawaiite", 
            "mafite", "mugearite", "oceanite", "palagonite", "tachylyte", "tholeiite", 
            "mafic", "melaphyre", "greenstone", "spilite", "greenschist", "blueschist", 
            "basite", "metabasite",)
        andesite = ("andesit", "boninite", "icelandite", "marianite", "adakite", 
            "propylite",)
        dacite = ("dacit", "santorinite", "ignimbrite",)
        rhyolite = ("rhyolit", "felsite", "liparite", "felsic", "silicic", "pumice", 
            "obsidian", "dellenite", "rhyodacite", "ignimbrite", "lenticulite", 
            "halleflinta", "leptite",)
        alk_volc = ( "polzenite", "hauynite", "arsoite", "benmoreite", "camptonite", 
            "ciminite", "damkjernite", "damtjernite", "domite", "fortunite", "gauteite",
            "kenyte", "keratophyre", "kersantite", "kivite", "lampro", "madupite",
            "minette", "monchiquite", "mondhaldeite", "orendite", "phonolite", "sannaite", 
            "trachyte", "wyomingite", "fonolito", "tinguaite",  "ordanchite", "melilit", 
            "katungite", "vsbergite", "analcimite", "ankaratrite", "etindite", "foidite", 
            "grazinite", "hauynophyre", "kalsilit", "leucitite", "mafurite", "melafoidite",
            "nephelinite","ugandite", "ottajanite", "melnoite", "pantellerite", "comendite", 
            "latite", "tristanite", "augitite", "absarokite", "shoshonite", "linosaite",
            "bergalite", "alnoite", "kimberlite",  "orangeite", "diatreme", "pipe", )
        volcaniclast = ("tonstein", "peperite", "volcaniclastic", "lahar",)
        volc = ("volcanic", "extrusive", "lava", "eutaxite", "vitrophyre", "volcan", 
            "ash", "ashfall", "tuff",  "tephra", "cinder", "porphyrite", 
            "vulsinite", "glass", "vitrophere",)
            
        # Plutonic
        peridotite = ("olivinite", "dunit", "lherzolite", "peridot", "harzburg", 
             "wehrlite", "serpentin", "soapstone", "talc", "alkremite", )
        pyroxenite = ("bronzitite", "pyroxenite", "enstatitite", "websterite",
            "hornblendite", "cortlandite",)
        gabbro = ("gabbro", "mafraite", "allivalite", "anorthosite", "diabase", "dolerit", 
            "leucophyre", "glenmuirite", "jotunite", "labradorite", "luscladite", 
            "theralite", "norite", "troctolite", "sebastianite", "eclogite", "amphibolit", 
            "rodingite", "corganite", "corgaspinite",)
        diorite = ("diorit", "jotunite", "marscoite", "sanukite",)
        trondhjemite =  ("trondhjemite", "trond",)
        tonalite = ("tonalit", "adamellite", "enderbite", "enderbite",)
        granodiorite = ("granodiorite")
        granite = ("granit", "microgranite", "adamellite", "aplite", "charnockite", 
            "granophyre", "rapakivi", "monzonit", "mangerite", "greisen", "pegmat",)
        alk_plut = ("syenit", "alaskite", "borolanite", "bostonite", "durbachite", 
            "foyaite", "jacupirangite", "juvite", "kentallenite", "larvikite", "lujavrite",
            "nordmarkite", "orthoclasite", "shonkinite", "sommaite", "kaersutitite",
            "lestiwarite", "puglianite", "vaugnerite","fergusite", "ijolite", 
            "melteigite", "missourite", "tannbuschite", "buchonite", "campanite", 
            "murambite", "tephrite", "tahitite", "vicoite", "urtite", "ankaramite", 
            "basanit", "limburgite", "biotitite", "riedenite", "glimmerite", "kamafugite",
            "turjaite", "essexite", "yamaskite", "teschenite", "crinanite", "vibetoite",  
            "uncompahgrite", "apatitite", "nelsonite", "phoscorite", "kullaite", )
        plut = ("plutonic", "pluton", "intrusive", "intrus", "sill", "dike", "stock", 
            "laccolith", "lopolith", "batholith", "porphyry", "megacryst",
            "hypabyssal", "chromitite", "topazite", )
            
        # Undefined igneous
        carbonatite = ("alvikite", "carbonatite", "beforsite", "rauhaugite", "sovite",
            "breunneritite",)
        ign = ("igneous", "metaign", "orthogneiss", "ortho", "meta-ign", "zeolite",)

        # Undefined metamorphic
        met = ("crystalline", "migma", "alter", "hydrothermal", "basement", 
            "high grade metamorphic", "meta", "granulit", "granofels", "gneiss", "schist", 
            "hornfels", "garnet", "buchite", "epidot", "fenite", "albitite", "chloritite", 
            "phlogopitite", "sericitite", "tactite", "tourmalinite", "unakite", 
            "vogesite", "gossan", "palagonite", "sanidinite",)

        # Cover
        cover = ("lluv", "fluv", "boulder", "gravel", "aleurite", "glaci", "till", "loess", 
            "regolith", "debris", "fill", "slide", "unconsolidated", "talus", "stream", 
            "beach", "terrace", "placer", "paleosol", "mass-wasting", "pebble", "cover", 
            "quaternary", "soil", "laterite", "surficial deposits", "scree", "peat", 
            "swamp", "marsh", "water", "ice")

        if inclusive || major
            # Major types include all minor subtypes
            sed = (siliciclast..., shale..., carb..., evap..., chert..., 
                phosphorite..., coal..., sed...,)
            ign = (komatiite..., basalt..., andesite..., dacite..., rhyolite..., 
                alk_volc..., volcaniclast..., volc..., peridotite..., pyroxenite..., 
                gabbro..., diorite..., trondhjemite..., tonalite..., granodiorite..., 
                granite..., alk_plut..., plut..., )

            # If only returning major types, do that 
            if major
                return (sed=sed, ign=ign, met=met, cover=cover)
            end
        end

        # If returning minor types, do that. The major types have been set to be inclusive
        # in the if statement above, if we wanted that
        return (
            # Sedimentary
            siliciclast=siliciclast, shale=shale, carb=carb, evap=evap, chert=chert, 
            phosphorite=phosphorite, coal=coal, sed=sed,

            # Volcanic
            komatiite=komatiite, basalt=basalt, andesite=andesite, dacite=dacite, 
            rhyolite=rhyolite, alk_volc=alk_volc, volcaniclast=volcaniclast, volc=volc, 
            
            # Plutonic
            peridotite=peridotite, pyroxenite=pyroxenite, gabbro=gabbro, diorite=diorite, 
            trondhjemite=trondhjemite, tonalite=tonalite, granodiorite=granodiorite, 
            granite=granite, alk_plut=alk_plut, plut=plut, 

            # Metamorphic
            met=met,

            # Cover
            cover=cover,
        )
    end


    """
    ```julia
    nondescriptive()
    ```

    Get a list of metamorphic rock names that do not provide useful information about the 
    geochemical composition of the rock.

    """
    function nondescriptive()
        return ("meta", "buchite", "tactite","greenschist", "alter", "hydrothermal",
            "crystalline", "basement", "skarn", "schist", "gneiss", "granulit", 
            "granofels", "sanidinite", "migma", "high grade metamorphic", "mylonit", 
            "cataclasite", "melange", "gouge", "tecton",
        )
    end

    """
    ```julia
    get_cats(major::Bool, npoints::Int64)
    ```
    
    Meow! 🐈 Initialize a NamedTuple of `npoints`-element BitVectors for each defined rock 
    type.
    
    If `major` is `true`, only major types are defined. If `major` is `false`, all subtypes
    are defined.

    See also: `get_minor_types` and `get_rock_class`.

    # Example
    ```julia
    typelist, cats = get_cats(true, 10)
    ```
    """
    function get_cats(major::Bool, npoints::Int64)
        typelist = get_rock_class(major=major)

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