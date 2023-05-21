# this is SIGNIFICANTLY worse than broadcasting...
function match_rocktype2(rocktype, rockname, rockdescrip; major=false)
    npoints = length(rocktype)

    # Rock types based on EarthChem system in RockNameInference.m
    # Sedimentary
    siliciclasttypes = ["siliciclast", "conglomerat", "sand", "psamm", "arenit", "arkos", "silt",]
    shaletypes = ["mud", "clay","shale", "wacke", "argillite", "argillaceous", "flysch", "pelit", 
        "turbidite",]
    carbtypes = ["carbonate", "limestone", "dolo", "marl", "chalk", "travertine", "tavertine", 
        "teravertine", "tufa",]
    cherttypes = ["chert", "banded iron",]
    evaporitetypes = ["evaporite", "gypsum", "salt flat",]
    coaltypes = ["coal", "anthracite",]
    
    sedtypes = vcat(["sediment", "fluv", "clast", "gravel", "pebble", "caliche", "boulder", 
        "diamict", "tillite", "stream", "beach", "terrace",  "marine deposits",  "paleosol"],
        siliciclasttypes, shaletypes, carbtypes, cherttypes, evaporitetypes, coaltypes)

    # Igneous
    volctypes = ["volcan", "lava", "lahar", "ignimbrite", "ashfall", "tuff", "diatreme",
        "pipe", "basalt", "andesit", "dacit", "rhyolit", "pillow", "carbonatite", "tephra", 
        "obsidian", "ash", "scoria", "pumice", "cinder", "latite", "basanite", "phonolite", 
        "fonolito", "trachyte", "palagonite", "mugearite", "kimberlite", "ultramafitite", 
        "komatiite",]
    hypabyssaltypes = ["intrus", "hypabyssal", "sill", "dike", "stock", "laccolith", "lopolith", 
        "dolerit", "diabase", "porphyry", "microgranite"]
    pluttypes = vcat(["pluton", "batholith", "granit", "tonalit", "gabbro", "norite", "diorit", 
        "monzonit", "syenit", "peridot", "dunit", "harzburg", "anorthosite", "mangerite", 
        "charnockite", "pegmatite", "aplite", "trond", "essexite", "pyroxenite", "adamellite", 
        "porphyry", "megacryst", "rapakivi", "bronzitite", "alaskite", "troctolite",], 
        hypabyssaltypes)
    
    igntypes = vcat(["igneous", "silicic ", "mafic", "felsic", "basite",], volctypes, pluttypes)

    # Metamorphic
    metasedtypes = ["para", "metased", "meta-sed", "quartzite", "marble", "slate", "phyllite",]
    metaigntypes = ["ortho", "metaign", "meta-ign", "serpentin", "amphibolit", "greenstone",
        "eclogite", "metabasite",]
    lowgradetypes = ["slate", "phyllite", "serpentin", "greenstone", "greenschist", "zeolite", 
        "gossan", "alter", "hydrothermal", "palagonite",]
    highgradetypes = ["crystalline", "basement", "marble", "skarn", "schist", "blueschist", 
        "gneiss", "amphibolit", "eclogite", "granulit", "hornfels", "granofels", "sanidinite", 
        "migma", "enderbite", "anorthosite", "charnockite", "pyroxenite", "peridot", "dunit", 
        "harzburg", "high grade metamorphic"]
    cataclastictypes = ["mylonit", "cataclasite", "melange", "gouge", "tecton",]
    
    mettypes = vcat(["meta", "calc silicate",], metasedtypes, metaigntypes, lowgradetypes, 
        highgradetypes, cataclastictypes)

    # Cover
    covertypes = ["cover", "unconsolidated", "quaternary", "lluv", "soil", "regolith", 
        "laterite", "surficial deposits", "talus", "scree", "mass-wasting", "slide", 
        "peat", "swamp", "marsh", "water", "ice", "glaci", "till", "loess", "gravel", 
        "debris"]
    
    # Define what we're looking for
    if major
        # Preallocate
        cats = (
            sed = falses(npoints),
            ign = falses(npoints),
            met = falses(npoints),
            cover = falses(npoints)
        )
        typelist = [sedtypes, igntypes, mettypes, covertypes]
    else
        cats = (
            siliciclast = falses(npoints),
            shale = falses(npoints),
            carb = falses(npoints),
            chert = falses(npoints),
            evaporite = falses(npoints),
            coal = falses(npoints),
            sed = falses(npoints),

            volc = falses(npoints),
            plut = falses(npoints),
            ign = falses(npoints),

            metased = falses(npoints),
            metaign = falses(npoints),
            met = falses(npoints),

            cover = falses(npoints),
        )
        typelist = [siliciclasttypes, shaletypes, carbtypes, cherttypes, evaporitetypes, 
            coaltypes, sedtypes, volctypes, pluttypes, igntypes, metasedtypes, metaigntypes, 
            mettypes, covertypes
        ]
    end

    # Check major lithology first
    for j in eachindex(typelist)
        for i = 1:length(typelist[j])
            for k in 1:length(cats[j])
                cats[j][k] |= match(r"major.*?{(.*?)}", rocktype[k]) |> x -> isa(x,RegexMatch) ? containsi(x[1], typelist[j][i]) : false
            end
        end
    end

    # Check the rest of rocktype
    not_matched = .!(cats.sed .| cats.ign .| cats.met .| cats.cover)
    # wtf in this uses 64 MB of memory... has to be >5MB to be competative
    for j in eachindex(typelist)
        for i = 1:length(typelist[j])
            for k in 1:length(cats[j])
                if cats[j][k]
                    continue
                else
                    cats[j][k] = containsi(rocktype[k], typelist[j][i])     # 656 bytes, vs broadcasting is 19.7 KB...
                end
            end
        end
    end

    # Then rockname
    not_matched = .!(cats.sed .| cats.ign .| cats.met .| cats.cover)
    for j in eachindex(typelist)
        for i = 1:length(typelist[j])
            for k in 1:length(cats[j][not_matched])
                cats[j][not_matched][k] |= containsi(rockname[not_matched][k], typelist[j][i])
            end
        end
    end

    # Then rockdescrip
    not_matched = .!(cats.sed .| cats.ign .| cats.met .| cats.cover)
    for j in eachindex(typelist)
        for i = 1:length(typelist[j])
            for k in 1:length(cats[j][not_matched])
                cats[j][not_matched][k] |= containsi(rockdescrip[not_matched][k], typelist[j][i])
            end
        end
    end

    return cats
end

# This could potentially be made more efficient with a loop...
# but even at 50000 points it still only uses about 6KB, so not a priority
function find_unmatched(cats; major=false)
    if major
        return .!(cats.sed .| cats.ign .| cats.met .| cats.cover)
    else
        cats.sed .= cats.sed .| cats.siliciclast .| cats.shale .| cats.carb .| cats.chert .| 
            cats.evaporite .| cats.coal
        cats.ign .= cats.ign .| cats.volc .| cats.plut
        cats.met .= cats.met .| cats.metased .| cats.metaign

        return .!(cats.sed .| cats.ign .| cats.met .| cats.cover)
    end
end

# relatively consistently slower, but with slightly fewer allocations
function match_reduced(rocktype, rockname, rockdescrip)
    npoints = length(rocktype)

    # Rock types based on EarthChem system in RockNameInference.m
        # Sedimentary
        siliciclasttypes = ["siliciclast", "conglomerat", "sand", "psamm", "arenit", "arkos", "silt",]
        shaletypes = ["mud", "clay","shale", "wacke", "argillite", "argillaceous", "flysch", "pelit", 
            "turbidite",]
        carbtypes = ["carbonate", "limestone", "dolo", "marl", "chalk", "travertine", "tavertine", 
            "teravertine", "tufa",]
        cherttypes = ["chert", "banded iron",]
        evaporitetypes = ["evaporite", "gypsum", "salt flat",]
        coaltypes = ["coal", "anthracite",]
        
        sedtypes = vcat(["sediment", "fluv", "clast", "gravel", "pebble", "caliche", "boulder", 
            "diamict", "tillite", "stream", "beach", "terrace",  "marine deposits",  "paleosol"],
            siliciclasttypes, shaletypes, carbtypes, cherttypes, evaporitetypes, coaltypes)

        # Igneous
        volctypes = ["volcan", "lava", "lahar", "ignimbrite", "ashfall", "tuff", "diatreme",
            "pipe", "basalt", "andesit", "dacit", "rhyolit", "pillow", "carbonatite", "tephra", 
            "obsidian", "ash", "scoria", "pumice", "cinder", "latite", "basanite", "phonolite", 
            "fonolito", "trachyte", "palagonite", "mugearite", "kimberlite", "ultramafitite", 
            "komatiite",]
        hypabyssaltypes = ["intrus", "hypabyssal", "sill", "dike", "stock", "laccolith", "lopolith", 
            "dolerit", "diabase", "porphyry", "microgranite"]
        pluttypes = vcat(["pluton", "batholith", "granit", "tonalit", "gabbro", "norite", "diorit", 
            "monzonit", "syenit", "peridot", "dunit", "harzburg", "anorthosite", "mangerite", 
            "charnockite", "pegmatite", "aplite", "trond", "essexite", "pyroxenite", "adamellite", 
            "porphyry", "megacryst", "rapakivi", "bronzitite", "alaskite", "troctolite",], 
            hypabyssaltypes)
        
        igntypes = vcat(["igneous", "silicic ", "mafic", "felsic", "basite",], volctypes, pluttypes)

        # Metamorphic
        metasedtypes = ["para", "metased", "meta-sed", "quartzite", "marble", "slate", "phyllite",]
        metaigntypes = ["ortho", "metaign", "meta-ign", "serpentin", "amphibolit", "greenstone",
            "eclogite", "metabasite",]
        lowgradetypes = ["slate", "phyllite", "serpentin", "greenstone", "greenschist", "zeolite", 
            "gossan", "alter", "hydrothermal", "palagonite",]
        highgradetypes = ["crystalline", "basement", "marble", "skarn", "schist", "blueschist", 
            "gneiss", "amphibolit", "eclogite", "granulit", "hornfels", "granofels", "sanidinite", 
            "migma", "enderbite", "anorthosite", "charnockite", "pyroxenite", "peridot", "dunit", 
            "harzburg", "high grade metamorphic"]
        cataclastictypes = ["mylonit", "cataclasite", "melange", "gouge", "tecton",]
        
        mettypes = vcat(["meta", "calc silicate",], metasedtypes, metaigntypes, lowgradetypes, 
            highgradetypes, cataclastictypes)

        # Cover
        covertypes = ["cover", "unconsolidated", "quaternary", "lluv", "soil", "regolith", 
            "laterite", "surficial deposits", "talus", "scree", "mass-wasting", "slide", 
            "peat", "swamp", "marsh", "water", "ice", "glaci", "till", "loess", "gravel", 
            "debris"]

    cats = (
        siliciclast = falses(npoints),
        shale = falses(npoints),
        carb = falses(npoints),
        chert = falses(npoints),
        evaporite = falses(npoints),
        coal = falses(npoints),
        sed = falses(npoints),

        volc = falses(npoints),
        plut = falses(npoints),
        ign = falses(npoints),

        metased = falses(npoints),
        metaign = falses(npoints),
        met = falses(npoints),

        cover = falses(npoints),
    )
    typelist = [siliciclasttypes, shaletypes, carbtypes, cherttypes, evaporitetypes, 
        coaltypes, sedtypes, volctypes, pluttypes, igntypes, metasedtypes, metaigntypes, 
        mettypes, covertypes
    ]

    for j in eachindex(typelist)
        for i = 1:length(typelist[j])
            for k in 1:length(cats[j])
                cats[j][k] |= match(r"major.*?{(.*?)}", rocktype[k]) |> x -> isa(x,RegexMatch) ? containsi(x[1], typelist[j][i]) : false
            end
        end
    end

    # Check the rest of rocktype
    # wtf in this uses so much memory...
    # this isn't picking up anything? 
    not_matched = find_unmatched(cats, major=false)
    for j in eachindex(typelist)
        for i = 1:length(typelist[j])
            for k in 1:length(cats[j])
                if not_matched[k]
                    cats[j][k] |= containsi(rocktype[k], typelist[j][i])
                end
            end
        end
    end

    # Then rockname
    not_matched = find_unmatched(cats, major=false)
    for j in eachindex(typelist)
        for i = 1:length(typelist[j])
            for k in 1:length(cats[j])
                if not_matched[k]
                    cats[j][k] |= containsi(rockname[k], typelist[j][i])
                end
            end
        end
    end

    # Then rockdescrip
    not_matched = find_unmatched(cats, major=false)
    for j in eachindex(typelist)
        for i = 1:length(typelist[j])
            for k in 1:length(cats[j])
                if not_matched[k]
                    cats[j][k] |= containsi(rockdescrip[k], typelist[j][i])
                end
            end
        end
    end

    return cats
end

function match_reduced2(rocktype, rockname, rockdescrip)
    npoints = length(rocktype)

    # Rock types based on EarthChem system in RockNameInference.m
        # Sedimentary
        siliciclasttypes = ["siliciclast", "conglomerat", "sand", "psamm", "arenit", "arkos", "silt",]
        shaletypes = ["mud", "clay","shale", "wacke", "argillite", "argillaceous", "flysch", "pelit", 
            "turbidite",]
        carbtypes = ["carbonate", "limestone", "dolo", "marl", "chalk", "travertine", "tavertine", 
            "teravertine", "tufa",]
        cherttypes = ["chert", "banded iron",]
        evaporitetypes = ["evaporite", "gypsum", "salt flat",]
        coaltypes = ["coal", "anthracite",]
        
        sedtypes = vcat(["sediment", "fluv", "clast", "gravel", "pebble", "caliche", "boulder", 
            "diamict", "tillite", "stream", "beach", "terrace",  "marine deposits",  "paleosol"],
            siliciclasttypes, shaletypes, carbtypes, cherttypes, evaporitetypes, coaltypes)

        # Igneous
        volctypes = ["volcan", "lava", "lahar", "ignimbrite", "ashfall", "tuff", "diatreme",
            "pipe", "basalt", "andesit", "dacit", "rhyolit", "pillow", "carbonatite", "tephra", 
            "obsidian", "ash", "scoria", "pumice", "cinder", "latite", "basanite", "phonolite", 
            "fonolito", "trachyte", "palagonite", "mugearite", "kimberlite", "ultramafitite", 
            "komatiite",]
        hypabyssaltypes = ["intrus", "hypabyssal", "sill", "dike", "stock", "laccolith", "lopolith", 
            "dolerit", "diabase", "porphyry", "microgranite"]
        pluttypes = vcat(["pluton", "batholith", "granit", "tonalit", "gabbro", "norite", "diorit", 
            "monzonit", "syenit", "peridot", "dunit", "harzburg", "anorthosite", "mangerite", 
            "charnockite", "pegmatite", "aplite", "trond", "essexite", "pyroxenite", "adamellite", 
            "porphyry", "megacryst", "rapakivi", "bronzitite", "alaskite", "troctolite",], 
            hypabyssaltypes)
        
        igntypes = vcat(["igneous", "silicic ", "mafic", "felsic", "basite",], volctypes, pluttypes)

        # Metamorphic
        metasedtypes = ["para", "metased", "meta-sed", "quartzite", "marble", "slate", "phyllite",]
        metaigntypes = ["ortho", "metaign", "meta-ign", "serpentin", "amphibolit", "greenstone",
            "eclogite", "metabasite",]
        lowgradetypes = ["slate", "phyllite", "serpentin", "greenstone", "greenschist", "zeolite", 
            "gossan", "alter", "hydrothermal", "palagonite",]
        highgradetypes = ["crystalline", "basement", "marble", "skarn", "schist", "blueschist", 
            "gneiss", "amphibolit", "eclogite", "granulit", "hornfels", "granofels", "sanidinite", 
            "migma", "enderbite", "anorthosite", "charnockite", "pyroxenite", "peridot", "dunit", 
            "harzburg", "high grade metamorphic"]
        cataclastictypes = ["mylonit", "cataclasite", "melange", "gouge", "tecton",]
        
        mettypes = vcat(["meta", "calc silicate",], metasedtypes, metaigntypes, lowgradetypes, 
            highgradetypes, cataclastictypes)

        # Cover
        covertypes = ["cover", "unconsolidated", "quaternary", "lluv", "soil", "regolith", 
            "laterite", "surficial deposits", "talus", "scree", "mass-wasting", "slide", 
            "peat", "swamp", "marsh", "water", "ice", "glaci", "till", "loess", "gravel", 
            "debris"]

    cats = (
        siliciclast = falses(npoints),
        shale = falses(npoints),
        carb = falses(npoints),
        chert = falses(npoints),
        evaporite = falses(npoints),
        coal = falses(npoints),
        sed = falses(npoints),

        volc = falses(npoints),
        plut = falses(npoints),
        ign = falses(npoints),

        metased = falses(npoints),
        metaign = falses(npoints),
        met = falses(npoints),

        cover = falses(npoints),
    )
    typelist = [siliciclasttypes, shaletypes, carbtypes, cherttypes, evaporitetypes, 
        coaltypes, sedtypes, volctypes, pluttypes, igntypes, metasedtypes, metaigntypes, 
        mettypes, covertypes
    ]

    # Check major lithology first
    for j in eachindex(typelist)
        for i = 1:length(typelist[j])
            cats[j] .|= (match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], typelist[j][i]) : false)
        end
    end

    # Check the rest of rocktype
    not_matched = find_unmatched(cats, major=false)
    for j in eachindex(typelist)
        for i = 1:length(typelist[j])
            cats[j][not_matched] .|= containsi.(rocktype[not_matched], typelist[j][i])
        end
    end

    # Then rockname
    not_matched = find_unmatched(cats, major=false)
    for j in eachindex(typelist)
        for i = 1:length(typelist[j])
            cats[j][not_matched].|= containsi.(rockname[not_matched], typelist[j][i])
        end
    end

    # Then rockdescrip
    not_matched = find_unmatched(cats, major=false)
    for j in eachindex(typelist)
        for i = 1:length(typelist[j])
            cats[j][not_matched].|= containsi.(rockdescrip[not_matched], typelist[j][i])
        end
    end

    return cats
end


## Test functions
    t = match_reduced(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)
    t2 = match_reduced2(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)

    @test t.siliciclast == t2.siliciclast
    @test t.shale == t2.shale
    @test t.carb == t2.carb
    @test t.chert == t2.chert
    @test t.evaporite == t2.evaporite
    @test t.coal == t2.coal
    @test t.sed == t2.sed

    @test t.volc == t2.volc
    @test t.plut == t2.plut
    @test t.ign == t2.ign

    @test t.metased == t2.metased
    @test t.metaign == t2.metaign
    @test t.met == t2.met

    @test t == t2