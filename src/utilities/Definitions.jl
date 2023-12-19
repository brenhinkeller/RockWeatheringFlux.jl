# Definitions, file names, etc. Contains most hard-coded values.

## --- File names

    # 500 Macrostrat samples
        # macrostrat_io = "output/N_500/toy_responses.h5"
        # matchedbulk_io = "output/N_500/toy_bulkidx.tsv"

        # eroded_out = "output/N_500/toy_erodedmaterial.h5"
        # ucc_out = "results/toy_exposedcrust.tsv"
        # erodedabs_out = "results/toy_erodedmaterial_abs.csv"
        # erodedrel_out = "results/toy_erodedmaterial_rel.csv"

    # 250_000 Macrostrat samples
        macrostrat_io = "output/N_250K/250K_responses.h5"
        matchedbulk_io = "output/N_250K/250K_bulkidx.tsv"

        eroded_out = "output/N_250K/250K_erodedmaterial.h5"
        ucc_out = "results/250K_exposedcrust.tsv"
        erodedabs_out = "results/250K_erodedmaterial_abs.tsv"
        erodedrel_out = "results/250K_erodedmaterial_rel.tsv"

    # 1_000_000 Macrostrat samples
        # macrostrat_io = "output/N_1M/1M_responses.h5"
        # matchedbulk_io = "output/N_1M/1M_bulkidx.tsv"

        # eroded_out = "output/N_1M/1M_erodedmaterial.h5"
        # ucc_out = "results/1M_exposedcrust.tsv"
        # erodedabs_out = "results/1M_erodedmaterial_abs.csv"
        # erodedrel_out = "results/1M_erodedmaterial_rel.csv"


## --- Color names for visualization

    # Gradients
    c_gradient = :jet1

    # Resampled
    c_rs = :grey

    # Observed samples by rock type
    lithclass = importdataset("data/lithclass_colors.tsv",'\t', importas=:Tuple)
    colortext = (
        ign = "Diabase", 
            volc = "Volcanic rock", 
                komatiite = "Ultramafitite",
                basalt = "Basalt",
                andesite = "Andesite",
                dacite = "Dacite",
                rhyolite = "Rhyolite",
                alk_volc = "Alkalic volcanic rock",
                volcaniclast = "Mixed volcanic/clastic rock",
            plut = "Plutonic rock",
                peridotite = "Peridotite",
                pyroxenite = "Pyroxenite",
                gabbro = "Gabbro",
                diorite = "Diorite",
                trondhjemite = "Trondhjemite",
                tonalite = "Tonalite",
                granodiorite = "Granodiorite",
                granite = "Granite",
                alk_plut = "Alkalic intrusive rock",
            carbonatite = "Intrusive carbonatite",
        sed = "Sedimentary rock", 
            siliciclast = "Medium-grained mixed clastic rock", 
            shale = "Shale", 
            carb = "Carbonate rock",
            chert = "Chert", 
            evap = "Evaporite", 
            coal = "Coal", 
            phosphorite = "Phosphorite",
        met = "Metamorphic rock",
        cover = "Unconsolidated material",
    )

    """
    Colors for each rock type from `get_rock_class`, based on the USGS lithologic 
    classification colors: https://mrdata.usgs.gov/catalog/lithclass-color.php.
    """
    colors = Dict{Symbol, RGB{Float64}}()
    for c in eachindex(colortext)
        i = findfirst(x -> x == colortext[c], lithclass.text)
        colors[c] = RGB(lithclass.r[i]/255, lithclass.g[i]/255, lithclass.b[i]/255)
    end
    colors = NamedTuple{Tuple(keys(colors))}(values(colors))

    
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

    """
    ```julia
    get_REEs()
    ```

    Rare earth elements. Note that Pm is often excluded from REE diagrams.
    """
    get_REEs() = return [:La, :Ce, :Pr, :Nd, :Pm, :Sm, :Eu, :Gd, :Tb, :Dy, :Ho, :Er, :Tm, 
        :Yb, :Lu,
    ]

    """
    ```julia
    get_chondrite_norm
    ```
    Chondrite values from Taylor and McLennan (1985).
    """
    get_chondrite_norm() = return (
        La = 0.367,
        Ce = 0.957,
        Pr = 0.137,
        Nd = 0.711,
        Sm = 0.231,
        Eu = 0.087,
        Gd = 0.306,
        Tb = 0.058,
        Dy = 0.381,
        Ho = 0.085,
        Er = 0.249,
        Tm = 0.036,
        Yb = 0.248,
        Lu = 0.038,
    )

    
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
            "breccia", "quartzite", "quarzite")
        shale = ("lutite", "mud", "clay", "shale", "wacke", "argillite", "argillaceous", 
            "flysch", "pelit", "turbidite", "tasmanite", "slate", "phyllite", 
            "metapellite", "micaschist", "mica schist")
        carb = ("carbonate", "calcite", "limestone", "dolo", "marl", "chalk", "coquina", 
            "biogenic", "travertine", "tavertine", "tufa", "calcarenite", "teravertine", 
            "marble", "calc silicate", "calcsilicate", "skarn",  "calcrete", "siderite", 
            "magnesite")
        evap = ("evaporite", "anhydrite", "gypsum", "trona", "halite", "sylvite", 
            "salt flat", "caliche", "exhalite", "sulfate deposit")
        chert = ("chert", "opal", "porcellanite", "diatomite", "novaculite", "iron", 
            "taconite", "banded iron", "spiculite",)
        phosphorite = ("phosphorite", "phosphate")
        coal = ("coal", "anthracite", "lignite", "bitumen")
        sed = ("sediment", "clast", "diamict", "tillite", "stream deposits", 
            "beach deposits", "terrace", "marine deposits", "paleosol", "glauconite", 
            "meta-sed", "metased", "paragneiss", "para", "melange")

        # Volcanic
        komatiite = ("komatiite", "meimechite", "ultramafitite")
        basalt = ("basalt", "pillow", "scoria", "picrite", "anamesite", "hawaiite", 
            "mafite", "mugearite", "oceanite", "palagonite", "tachylyte", "tholeiite", 
            "mafic", "melaphyre", "greenstone", "spilite", "greenschist", "blueschist", 
            "basite", "metabasite", "hyaloclastite", "tholleiite")
        andesite = ("andesit", "andeste", "boninite", "icelandite", "marianite", "adakite", 
            "propylite",)
        dacite = ("dacit", "santorinite", "ignimbrite",)
        rhyolite = ("rhyolit", "felsite", "liparite", "felsic", "silicic", "pumice", 
            "obsidian", "dellenite", "rhyodacite", "ignimbrite", "lenticulite", 
            "halleflinta", "leptite", "rhyoite")
        alk_volc = ("polzenite", "hauynite", "arsoite", "benmoreite", "camptonite", 
            "ciminite", "damkjernite", "damtjernite", "dankjernite", "domite", "fortunite", 
            "gauteite","kenyte", "keratophyre", "kersantite", "kivite", "lampro", "madupite",
            "minette", "monchiquite", "mondhaldeite", "orendite", "phonolite", "sannaite", 
            "trachyte", "wyomingite", "fonolito", "tinguaite",  "ordanchite", "melilit", 
            "katungite", "vsbergite", "analcimite", "ankaratrite", "etindite", "foidite", 
            "grazinite", "hauynophyre", "kalsilit", "leucitite", "mafurite", "melafoidite",
            "nephelinite","ugandite", "ottajanite", "melnoite", "pantellerite", "comendite", 
            "latite", "tristanite", "augitite", "absarokite", "shoshonite", "linosaite",
            "bergalite", "alnoite", "alnã\u96ite", "kimberlite",  "orangeite", "diatreme", 
            "pipe", "alkaline volcan", "alkalic igneous", "malignite",)
        volcaniclast = ("tonstein", "peperite", "volcaniclastic", "lahar",)
        volc = ("volcanic", "extrusive", "lava", "eutaxite", "vitrophyre", "volcan", 
            "ash", "ashfall", "tuff",  "tephra", "cinder", "porphyrite", 
            "vulsinite", "glass", "vitrophere", "pyroclastic", "agglomerate")
            
        # Plutonic
        peridotite = ("periodotite", "olivinite", "dunit", "lherzolite", "peridot", "harzburg", 
             "wehrlite", "wehlerite", "serpentin", "soapstone", "talc", "alkremite",)
        pyroxenite = ("bronzitite", "pyroxenite", "enstatitite", "websterite",
            "hornblendite", "cortlandite",)
        gabbro = ("gabbro", "gabro", "mafraite", "allivalite", "anorthosite", "diabase", 
            "dolerit", "leucophyre", "glenmuirite", "jotunite", "labradorite", "luscladite", 
            "theralite", "norite", "troctolite", "sebastianite", "eclogite", "amphibolit", 
            "rodingite", "corganite", "corgaspinite",)
        diorite = ("diorit", "iorite", "jotunite", "marscoite", "sanukite",)
        trondhjemite =  ("trondhjemite", "trond",)
        tonalite = ("tonalit", "adamellite", "enderbite", "enderbite",)
        granodiorite = ("granodiorite",)
        granite = ("granit", "microgranite", "adamellite", "aplite", "charnockite", 
            "granophyre", "rapakivi", "monzonit", "monzonize", "mangerite", "greisen", 
            "pegmat", "adamelllite", "adamelite", "unakite")
        alk_plut = ("syenit", "seyenite", "alaskite", "borolanite", "bostonite", "durbachite", 
            "foyaite", "jacupirangite", "juvite", "kentallenite", "larvikite", "lujavrite",
            "nordmarkite", "orthoclasite", "shonkinite", "sommaite", "kaersutitite",
            "lestiwarite", "puglianite", "vaugnerite","fergusite", "ijolite", "ljolite",
            "melteigite", "missourite", "tannbuschite", "buchonite", "campanite", 
            "murambite", "tephrite", "tahitite", "vicoite", "urtite", "ankaramite", 
            "basanit", "limburgite", "biotitite", "riedenite", "glimmerite", "kamafugite",
            "turjaite", "essexite", "yamaskite", "teschenite", "crinanite", "vibetoite",  
            "uncompahgrite", "apatitite", "nelsonite", "phoscorite", "kullaite", "malignite",
            "alkaline pluton")
        plut = ("plutonic", "pluton", "intrusive", "intrus", "sill", "dike", "stock", 
            "laccolith", "lopolith", "batholith", "porphyry", "megacryst",
            "hypabyssal", "chromitite", "topazite", )
            
        # Undefined igneous
        carbonatite = ("alvikite", "carbonatite", "beforsite", "rauhaugite", "sovite",
            "breunneritite", "fenite")
        ign = ("igneous", "metaign", "orthogneiss", "ortho", "meta-ign", "zeolite", "xenolith")

        # Undefined metamorphic
        met = ("crystalline", "migma", "alter", "hydrothermal", "basement", 
            "high grade metamorphic", "meta", "granulit", "granofels", "schist", "schsit", 
            "gneiss", "hornfels", "garnet", "spessartite", "melanite", "buchite", "epidot", 
            "fenite", "albitite", "chloritite", "phlogopitite", "sericitite", "tactite", 
            "tourmalinite", "vogesite", "gossan", "palagonite", "sanidinite", "mylonite")

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
        minorign = (:volc, :plut, :carbonatite)

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
                ign=ign,

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

    """
    ```julia
    rm_false_positives!(cats)
    ```

    Remove false matches from `cats`, where a false match is when a rock name used to
    identify samples is present in another rock name. For example, all grano*diorite* 
    samples will match with *diorite*.

    Currently removes:
     *  Diorite from grano*diorite*.
     *  Lignite from ma*lignite*.
    """
    function rm_false_positives!(cats)
        cats.diorite .&= cats.granodiorite  # Diorite / granodiorite
        cats.coal .&= .!cats.alk_volc       # Lignite / malignite 

        return cats
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