########################
## SubClonalSelection ##
########################

# state of the package 2018-08-01
# https://marcjwilliams1.github.io/quantifying-selection
# Pkg.clone("https://github.com/marcjwilliams1/SubClonalSelection.jl")

# Pkg.add("DataFrames", v"0.11.3")
# Pkg.pin("DataFrames", v"0.11.3")


# for parallel:
# julia -p 26

@everywhere using SubClonalSelection
@everywhere using DataFrames
@everywhere using CSV

srand(3458)

donors = ["joxm", "garx", "wahn", "vass", "ualf", "euts", "laey", "pipw", "oilg", "heja",
          "sehl", "feec", "gesg", "fikt", "vuna", "qonc", "xugn", "qolg", "puie", "fawm",
          "oaaz", "naju", "ieki", "rozh", "wetu", "nusw", "zoxy", "hipn", "lexy", "vils",
          "qayj", "kuco"]

filenames = readdir("data/AF")

subset = []
for donor in donors
    subset = vcat(subset, filter(x->contains(x, donor),filenames))
end

subset_petr = filter(x->contains(x, "petr"), subset)

@everywhere function runABCfitSubclone1(fileid)
    # run ABC estimation with parameters for Petr call set
    # detect max of 1 subclone
    mut_freq = CSV.read("data/AF/$fileid", delim="\t", types=Dict("chr"=>String))
    mut_freq[:af] = convert(Array{Float64}, mut_freq[:af])
    out = fitABCmodels(mut_freq[:af],
            splitext(basename(fileid))[1],
            resultsdirectory = "data/subclonal-output-1/",
            save = true,
            read_depth = 254,
            minreads = 3,
            fmin = 0.05,
            fmax = 0.45,
            maxiterations = 6*10^6,
            maxclones = 1, # only identify up to 1 subclone
            nparticles = 500,
            verbose = false,
            adaptpriors = false,
            ploidy = 2,
            mincellularity = 0.98)
end



@parallel for fileid = subset_petr
    print(fileid)
    runABCfitSubclone1(fileid)
end
