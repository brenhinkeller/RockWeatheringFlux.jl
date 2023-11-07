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
            chert = :deepskyblue, evap = :lightskyblue, coal = :midnightblue, 
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
    get_rock_class([major::Bool])
    ```

    Define sedimentary, igneous, and metamorphic rock types and subtypes. Metamorphic rocks
    are grouped with their protoliths when possible. Return a list of the minor types which
    map to each major type. 

    Define `major=true` to return only sedimentary, igneous, and metamorphic types. Rather
    than a list of subtypes, this will return an empty set of values.

    ### Possible subtypes
     *  Sedimentary: siliciclastic, shale, carbonate, evaporite, chert, phosphorite, coal,
        sedimentary (uncategorized)
     *  Igneous: volcanic, plutonic, carbonatite, igneous (uncategorized)
         *  Volcanic: komatiite, basalt, andesite, dacite, rhyolite, alkaline, 
            volcaniclastic, volcanic (uncategorized)
         *  Plutonic: peridotite, pyroxenite, gabbro, diorite, trondhjemite, tonalite, 
            granodiorite, granite, alkaline, plutonic (uncategorized)
      * Metamorphic: metamorphic (uncategorized)

    Subtypes are defined as a list of each bullet point, but the major (uncategorized) type
    is not included. That is, `minorsed` includes `siliciclastic`, `shale`, `carbonate`, 
    `evaporite`, `chert`, `phosphorite`, and `coal`, but not `sed`.

    The `minorvolc` and `minorplut` lists do not include `volc` or `plut`, but these are 
    included as subtypes of `ign`.

    ### Subtype Inclusion
    The lists of rock names matched to major rock types (e.g., sedimentary, igneous, etc.)
    do __not__ include the rock names which match to minor types. That is, the list of 
    sedimentary rock names does not include names which describe shales or cherts.

    Volcanic and plutonic rocks are similar. The list of volcanic rock names does not
    include rock names which describe basalts or andesites.

    To include these subtypes in the list of all rock names, set `major=true`.

    # Examples
    ```julia-repl
    julia> typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();

    julia> typelist, = get_rock_class(major=true);
    ```

    """
    function get_rock_class(; major::Bool=false)
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
        alk_volc = ("polzenite", "hauynite", "arsoite", "benmoreite", "camptonite", 
            "ciminite", "damkjernite", "damtjernite", "domite", "fortunite", "gauteite",
            "kenyte", "keratophyre", "kersantite", "kivite", "lampro", "madupite",
            "minette", "monchiquite", "mondhaldeite", "orendite", "phonolite", "sannaite", 
            "trachyte", "wyomingite", "fonolito", "tinguaite",  "ordanchite", "melilit", 
            "katungite", "vsbergite", "analcimite", "ankaratrite", "etindite", "foidite", 
            "grazinite", "hauynophyre", "kalsilit", "leucitite", "mafurite", "melafoidite",
            "nephelinite","ugandite", "ottajanite", "melnoite", "pantellerite", "comendite", 
            "latite", "tristanite", "augitite", "absarokite", "shoshonite", "linosaite",
            "bergalite", "alnoite", "kimberlite",  "orangeite", "diatreme", "pipe", 
            "alkaline volcan")
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
        granodiorite = ("granodiorite",)
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
            "uncompahgrite", "apatitite", "nelsonite", "phoscorite", "kullaite", 
            "alkaline pluton")
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

        # Define major and minor rock types
        minorsed = (:siliciclast, :shale, :carb, :evap, :chert, :phosphorite, :coal,)
        minorvolc = (:komatiite, :basalt, :andesite, :dacite, :rhyolite, :alk_volc, 
            :volcaniclast,)
        minorplut = (:peridotite, :pyroxenite, :gabbro, :diorite, :trondhjemite, :tonalite, 
            :tonalite, :granodiorite, :granite, :alk_plut,)
        minorign = (:volc, plut, :carbonatite)

        if major
            typelist = (sed=(minorsed, sed...,), ign=(minorign, ign...,), met=met, cover=cover)
            minors = ()
        else
            typelist = (
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

                # Igneous
                carbonatite=carbonatite,
                ign=carbonatite,

                # Metamorphic
                met=met,

                # Cover
                cover=cover,
            )
            minors = (minorsed, minorvolc, minorplut, minorign)
        end

        return typelist, minors...
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
        typelist, = get_rock_class(major=major)

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
        @warn "get_minor_types has been deprecated. Use get_rock_class instead."

        types = get_cats(false, 1)[2]
        allkeys = collect(keys(types))

        sed = findfirst(==(:sed), allkeys)
        ign = findfirst(==(:ign), allkeys)
        met = findfirst(==(:met), allkeys)

        return Tuple(allkeys[1:sed-1]), Tuple(allkeys[sed+1:ign-1]), Tuple(allkeys[ign+1:met-1])
    end


## --- End of file