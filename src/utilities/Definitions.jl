## --- Dataset switches 

    # Geochemical 
    # dataset = "bulk"
    dataset = "gard"

    # Lithologic
    # N = "500"
    N = "250K"
    # N = "1M"
    

## --- File names

    # Geochemical data 
    geochem_fid = "output/" * dataset * ".h5"

    # Lithologic data
    macrostrat_io = "output/N_" * N * "/" * N * "_responses.h5"
    
    # Intermediate files 
    matchedbulk_io = "output/N_" * N * "/" * N * "_bulkidx_" * dataset * ".tsv"
    eroded_out = "output/N_250K/250K_erodedmaterial_" * dataset * ".h5"

    # Bulk continental crust
    ucc_out = "results/" * N * "_exposedcrust_" * dataset * ".tsv"
    # TO DO: error output

    # Eroded material 
    erodedabs_out = "results/" * N * "_eroded_absolute_" * dataset * ".tsv"
    erodedrel_out = "results/" * N * "_eroded_fraction_" * dataset * ".tsv"
    erodedcomp_out = "results/" * N * "_eroded_composition_" * dataset * ".tsv"

    erodedabs_out_err = "results/" * N * "_eroded_absolute_err_" * dataset * ".tsv"
    erodedrel_out_err = "results/" * N * "_eroded_fraction_err_" * dataset * ".tsv"
    erodedcomp_out_err = "results/" * N * "_eroded_composition_err_" * dataset * ".tsv"

    export geochem_fid, macrostrat_io
    export matchedbulk_io, eroded_out
    export ucc_out
    export erodedabs_out, erodedrel_out, erodedcomp_out
    export erodedabs_out_err, erodedrel_out_err, erodedcomp_out_err


## --- Color names for visualization

    # Gradients
    # c_gradient = :jet1
    c_gradient = :nipy_spectral
    export c_gradient

    # Resampled
    c_rs = :grey
    export c_rs

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
    export colors

    """
    ```julia
    display_colors()
    ```

    Display all colors in `colors` to the plots window.
    """
    function display_colors()
        typelist, =  get_rock_class()
        k = keys(typelist)
        for i in eachindex(k)
            t[i] = colors[k[i]]
        end
        display(t)
    end
    export display_colors

    
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
    (https://www.ngu.no/filearchive/NGUPublikasjoner/NGUnr\\_322\\_Bulletin\\_35\\_Faye\\_35\\_53.pdf).

    Note that Cr₂O₃, MnO, NiO, and P₂O₅ are minor elements calculated as wt.% element 
    oxide.

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
    export get_elements


    """
    ```julia
    get_REEs()
    ```

    Rare earth elements. Note that Pm is often excluded from REE diagrams.
    """
    get_REEs() = return [:La, :Ce, :Pr, :Nd, :Pm, :Sm, :Eu, :Gd, :Tb, :Dy, :Ho, :Er, :Tm, 
        :Yb, :Lu,
    ]
    export get_REEs


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
            "breccia", "quartzite", "quarzite", "sst", "cgl")
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
            :granodiorite, :granite, :alk_plut,)
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
    export get_rock_class

    
    """
    ```julia
    get_rock_names()

    Major and minor rock names from `get_rock_class`, unabbreviated for plot labels.
    ```
    """
    get_rock_names() = return ("Siliciclastic", "Shale", "Carbonate", "Evaporite", 
        "Chert", "Phosphorite", "Coal", "Sedimentary", "Komatiite", 
        "Basalt", "Andesite", "Dacite", "Rhyolite", "Alkaline Volcanic", 
        "Volcaniclastic", "Volcanic", "Peridotite", "Pyroxenite", "Gabbro", 
        "Diorite", "Trondhjemite", "Tonalite", "Granodiorite", "Granite", 
        "Alkaline Plutonic", "Plutonic", "Carbonatite", "Igneous", 
        "Unspecified Metamorphic", "Cover"
    )
    export get_rock_names

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
    export get_cats


## --- End of file