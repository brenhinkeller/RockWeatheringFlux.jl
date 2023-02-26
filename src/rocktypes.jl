## --- Define rock types based on names or partial names of different rock types from GetBurwellBulkAge.m
  # Sedimentary
  sedtypes = ["sediment", "fluv", " clast", "siliciclast", "conglomerat", "gravel", "sand", "psamm", "arenit", "arkos", "silt", 
    "mud", "marl", "clay", "shale", "wacke", "argillite", "argillaceous", "pelit", "pebble", "carbonate", "limestone", "dolo", 
    "caliche", "chalk", "travertine", "tavertine", "teravertine", "tufa", "evaporite", " salt", "salt flat", "gypsum", "boulder", 
    "diamict", "tillite", "stream", "beach", "terrace", "chert", "banded iron", "coal", "anthracite", "marine deposits", "turbidite", 
    "flysch", "paleosol"]

  # Igneous
  volctypes = ["volcan", "lava", "lahar", "ignimbrite", "ashfall", "tuff", "diatreme", "pipe", "basalt", "andesit", "dacit", "rhyolit", 
    "pillow", "carbonatite", "tephra", "obsidian", "ash", "scoria", "pumice", "cinder", "latite", "basanite", "phonolite", "fonolito", 
    "trachyte", "palagonite", "mugearite", "kimberlite", "ultramafitite", "komatiite",]
  pluttypes = ["pluton", "batholith", "granit", "tonalit", "gabbro", "norite", "diorit", "monzonit", "syenit", "peridot", "dunit", 
    "harzburg", "anorthosite", "mangerite", "charnockite", "pegmatite", "aplite", "trond", "essexite", "pyroxenite", "adamellite", 
    "porphyry", "megacryst", "rapakivi", "bronzitite", "alaskite", "troctolite",]
  hypabyssaltypes = ["intrus", "hypabyssal", "sill", "dike", "stock", "laccolith", "lopolith", "dolerit", "diabase", "porphyry", 
    "microgranite"]
  igntypes = vcat(["igneous", "silicic ", "mafic", "felsic", "basite",],volctypes,pluttypes,hypabyssaltypes)

  # Metamorphic
  metasedtypes = ["para", "metased", "meta-sed", "schist", "quartzite", "marble", "skarn", "slate", "phyllite",]
  metaigntypes = ["ortho", "metaign", "meta-ign", "serpentin", "amphibolit", "greenstone", "eclogite", "metabasite", "migma",]
  mettypes = vcat(metasedtypes, metaigntypes, ["gneiss", "granulit", "hornfels", "granofels", "mylonit", "meta", "cataclasite", 
    "melange", "gouge", "tecton", "calc silicate"])
  lowgradetypes = ["slate", "phyllite", "serpentin", "greenstone", "greenschist", "zeolite", "gossan", "alter", "hydrothermal", 
    "palagonite",]
  highgradetypes = ["crystalline", "basement", "marble", "skarn", "blueschist", "gneiss", "amphibolit", "eclogite", "granulit", 
    "hornfels", "granofels", "sanidinite", "migma", "enderbite", "anorthosite", "charnockite", "pyroxenite", "peridot", "dunit", 
    "harzburg", "high grade metamorphic"]

  # Other
  covertypes = ["cover", "unconsolidated", "quaternary", "lluv", "soil", "regolith", "laterite", "surficial deposits", "talus", 
      "scree", "mass-wasting", "slide", "peat", "swamp", "marsh", "water", "ice", "glaci", "till", "loess", "gravel", "debris"]
  cataclastictypes = ["mylonit", "cataclasite", "melange", "gouge", "tecton",]