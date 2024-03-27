## --- Dataset switches 

    # Geochemical 
    # dataset = "bulk"
    # dataset = "gard"
    dataset = "combined"

    # Lithologic
    # N = "500"
    # N = "250K"
    N = "1M"
    
    if dataset == "bulk" || dataset == "gard"
        @error "EarthChem and Gard datasets may not be compatible with current code base."
    end

## --- File names

    # Geochemical data 
    geochem_fid = "output/" * dataset * ".h5"

    # Lithologic data
    macrostrat_io = "output/N_" * N * "/" * N * "_responses.h5"
    surficial_abundance_total_out = "results/" * N * "_surficial_abundance_mapped.tsv"
    surficial_abundance_out = "results/" * N * "_surficial_abundance.tsv"
    
    # Intermediate files 
    matchedbulk_io = "output/N_" * N * "/" * N * "_bulkidx_" * dataset * ".tsv"
    eroded_out = "output/N_250K/250K_erodedmaterial_" * dataset * ".h5"

    # Bulk continental crust
    ucc_out = "results/" * N * "_exposedcrust_" * dataset * ".tsv"
    ucc_out_err = "results/" * N * "_exposedcrust_err_" * dataset * ".tsv"

    # Eroded material: absolute, fractional contribution, composition
    erodedabs_out = "results/" * N * "_eroded_absolute_" * dataset * ".tsv"
    erodedrel_out = "results/" * N * "_eroded_fraction_" * dataset * ".tsv"
    erodedcomp_out = "results/" * N * "_eroded_composition_" * dataset * ".tsv"

    erodedabs_out_err = "results/" * N * "_eroded_absolute_err_" * dataset * ".tsv"
    erodedrel_out_err = "results/" * N * "_eroded_fraction_err_" * dataset * ".tsv"
    erodedcomp_out_err = "results/" * N * "_eroded_composition_err_" * dataset * ".tsv"

    export geochem_fid, macrostrat_io
    export surficial_abundance_total_out, surficial_abundance_out
    export matchedbulk_io, eroded_out
    export ucc_out, ucc_out_err
    export erodedabs_out, erodedrel_out, erodedcomp_out
    export erodedabs_out_err, erodedrel_out_err, erodedcomp_out_err    


## --- Color names

    # Gradients
    colorgradient = :nipy_spectral
    export colorgradient

    # Unified color palette 
    # colorpalette = [parse(Colorant, "#C84630"), parse(Colorant, "#286886"), 
    #     parse(Colorant, "#3E000C"), parse(Colorant, "#F5B700"), parse(Colorant, "#63A46C")
    # ]
    colorpalette = :berlin
    export colorpalette

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
        t = [colors[k] for k in keys(typelist)]
        display(t)
    end
    export display_colors

    
## --- Lists of elements: major, minor, REEs

    """
    ```julia
    get_elements()
    ```

    Define major and minor elements, returned as `Vector{Symbol}`s

    Major elements:
    * SiO₂, Al₂O₃, FeOT, TiO₂, MgO, CaO, Na₂O, K₂O, Volatiles

        
    Minor elements:
    * Ag, As, Au, B, Ba, Be, Bi, C, Cd, Ce, Cl, Co, Cr, Cs, Cu, Dy, Er, Eu, 
        F, Ga, Gd, Hf, Hg, Ho, I, In, Ir, La, Li, Lu, MnO, Mo, Nb, Nd, Ni, Os, P₂O₅, Pb, 
        Pd, Pt, Pr, Re, Rb, Sb, Sc, Se, S, Sm, Sn, Sr, Ta, Tb, Te, Th, Tl, Tm, U, V, W, Y, 
        Yb, Zn, Zr

    Major elements are in part defined based on Faye and Ødegård 1975 
    (https://www.ngu.no/filearchive/NGUPublikasjoner/NGUnr\\_322\\_Bulletin\\_35\\_Faye\\_35\\_53.pdf).

    Note that MnO and P₂O₅ are minor elements calculated as wt.% element 
    oxide.

    See also: `major_elements`

    """
    function get_elements()
        majors = [:SiO2,:Al2O3,:FeOT,:TiO2,:MgO,:CaO,:Na2O,:K2O,:Volatiles]
        minors = [:Ag,:As,:Au,:B,:Ba,:Be,:Bi,:C,:Cd,:Ce,:Cl,:Co,:Cr,:Cs,:Cu,
            :Dy,:Er,:Eu,:F,:Ga,:Gd,:Hf,:Hg,:Ho,:I,:In,:Ir,:La,:Li,:Lu,:MnO,:Mo,:Nb,:Nd,
            :Ni,:Os,:P2O5,:Pb,:Pd,:Pt,:Pr,:Re,:Rb,:Sb,:Sc,:Se,:S,:Sm,:Sn,:Sr,:Ta,:Tb,
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


## --- Litholgic class wordbanks

    """
    ```julia
    get_rock_class([major::Bool])
    ```

    Define sedimentary, igneous, and metamorphic rock classes and subclasses. Metamorphic rocks
    are grouped with their protoliths when possible. Return a list of the minor classes which
    map to each major class. 

    Define `major=true` to return only sedimentary, igneous, and metamorphic types. Rather
    than a list of subclasses, this will return an empty set of values.

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
            "magnesite", "boundstone", "packstone", "grainstone", )
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
        rhyolite = ("rhyolit", "felsite", "liparite", "silicic", "pumice", 
            "obsidian", "dellenite", "rhyodacite", "ignimbrite", "lenticulite", 
            "halleflinta", "rhyoite")
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
        trondhjemite =  ("trondhjemite", "trond", "ttg",)
        tonalite = ("tonalit", "adamellite", "enderbite", "ttg",)
        granodiorite = ("granodiorite", "ttg",)
        granite = ("granit", "microgranite", "adamellite", "aplite", "charnockite", 
            "granophyre", "rapakivi", "monzonit", "monzonize", "mangerite", "greisen", 
            "pegmat", "adamelllite", "adamelite", "unakite", "felsic",)
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
        met = ("crystalline", "migma", "alter", "hydrothermal", "basement", "leptite",
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
    get_metamorphic_class()
    ```

    As `get_rock_class`, but defines metasedimentary, metaigneous, and undifferentiated 
    metamorphic rock classes and subclasses. Does *not* return minor subclasses: see 
    `get_rock_class`.

    To match the format of `get_rock_class`, this function will return BitVectors for all
    defined minor classes, even if there are no rock names which match to metamorphic 
    rocks. Classes which will always be empty are:
      * `evaporite, phosphorite, coal, dacite, alkaline volcanic, volcaniclastics, 
        pyroxenite, diorite, trondhjemite, granodiorite, alk_plut, carbonatite, cover`

    """
    function get_metamorphic_class()
        # Sedimentary
        siliciclast = ("quartzite", "quarzite",)
        shale = ("pelit",  "slate", "phyllite", "metapellite", "micaschist", "mica schist")
        carb = ("marble", "calc silicate", "calcsilicate", "skarn")
        evap = ()
        chert = ("porcellanite",)
        phosphorite = ()
        coal = ()
        sed = ("meta-sed", "metased", "paragneiss", "para", "melange")

        # Volcanic
        komatiite = ("ultramafitite",)
        basalt = ("palagonite","greenstone", "spilite", "greenschist", "blueschist", 
            "basite", "metabasite", "melaphyre",)
        andesite = ("propylite",)
        dacite = ()
        rhyolite = ("halleflinta",)
        alk_volc = ()
        volcaniclast = ()
        volc = ()
            
        # Plutonic
        peridotite = ("serpentin", "soapstone", "talc", "alkremite")
        pyroxenite = ()
        gabbro = ("eclogite", "amphibolit", "leucophyre", "rodingite")
        diorite = ()
        trondhjemite =  ()
        tonalite = ("enderbite",)
        granodiorite = ()
        granite = ("greisen", "charnockite", "unakite")
        alk_plut = ()
        plut = ()
            
        # Undefined igneous
        carbonatite = ()
        ign = ("metaign", "orthogneiss", "ortho", "meta-ign", "zeolite",)

        # Undefined metamorphic
        met = ("crystalline", "migma", "alter", "hydrothermal", "basement", "leptite",
            "high grade metamorphic", "meta", "granulit", "granofels", "schist", "schsit", 
            "gneiss", "hornfels", "garnet", "spessartite", "melanite", "buchite", "epidot", 
            "fenite", "albitite", "chloritite", "phlogopitite", "sericitite", "tactite", 
            "tourmalinite", "vogesite", "gossan", "palagonite", "sanidinite", "mylonite")

        # Cover
        cover = ()

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

            # Igneous
            carbonatite=carbonatite,
            ign=ign,

            # Metamorphic
            met=met,

            # Cover
            cover=cover,
        )
    end
    export get_metamorphic_class
    
## --- Lithologic class accessory functions 

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

    """
    ```julia
    get_lithologic_class([matchedbulk_io], [macrostrat_io])
    ```

    Return standardized filters for lithologic class. Optionally specify paths for matched 
    samples and the Macrostrat responses file.

    # Return Values 
    * `match_cats`: lithologic class of samples, from the classes assigned during lithologic
        and geochemical sample matching. The `met` key is true if sample was mapped as a 
        metamorphic sample in Macrostrat, but no samples are exclusively metamorphic.
    * `metamorphic_cats`: true if the sample is mapped as metamorphic in Macrostrat.
    * `class`: as `match_cats`, but with an additional `bulk` key that is true for all samples.
    * `megaclass`: as `class`, but with additional filters for metasedimentary, metaigneous,
        and undifferentiated metamorphic samples.

    All filters are inclusive of minor lithologic classes.

    # Example
    ```julia
    match_cats, metamorphic_cats, class, megaclass = get_lithologic_class();
    ```

    """
    function get_lithologic_class(matchedbulk_io=matchedbulk_io, macrostrat_io=macrostrat_io)
        # Matched lithologic class
        fid = readdlm(matchedbulk_io)
        bulkidx = Int.(vec(fid[:,1]))
        t = @. bulkidx != 0
        match_cats = match_rocktype(string.(vec(fid[:,2]))[t]);

        # Metamorphic tag
        fid = h5open(macrostrat_io, "r")
        header = Tuple(Symbol.(read(fid["type"]["macro_cats_head"])))
        data = read(fid["type"]["metamorphic_cats"])
        data = @. data > 0
        metamorphic_cats = NamedTuple{header}([data[:,i][t] for i in eachindex(header)])
        close(fid)

        # Delete cover
        match_cats = delete_cover(match_cats)
        metamorphic_cats = delete_cover(metamorphic_cats)
        
        # Include minor types 
        include_minor!(match_cats)
        include_minor!(metamorphic_cats)

        # Set matched lithologic classes to include all identified metamorphic samples
        match_cats.met .|= (metamorphic_cats.sed .| metamorphic_cats.ign .| metamorphic_cats.met)

        # Allow indexing into all samples 
        class = merge(match_cats, (bulk=trues(length(match_cats[1])),))

        # As above, but with more metamorphic options
        metamorphic_cats.met .&= .!(metamorphic_cats.sed .| metamorphic_cats.ign)
        megaclass = merge(match_cats, (
            metased = metamorphic_cats.sed,
            metaign = metamorphic_cats.ign,
            met_undiff = metamorphic_cats.met,
            bulk=trues(length(match_cats[1])),
        ))

        return match_cats, metamorphic_cats, class, megaclass
    end
    export get_lithologic_class


## --- Matched sample metadata 
    
    """
    ```julia
    unique_sample(sample_ID, [p])
    ```

    Get the number of geochemical samples that explain `p`% of the matches. Default is 90%.

    """
    function unique_sample(sample_ID, p::Int=90)
        c = countmap(sample_ID)
        return count(<(percentile(values(c), p)), values(c))
    end
    export unique_sample


## --- End of file