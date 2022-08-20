module DrWatsonHelpers

import DrWatson: scriptsdir, datadir

export slurmscriptsdir, agatedatadir

slurmscriptsdir(args...) = scriptsdir("Slurm", args...)

agatedatadir(args...) = datadir("Agate", args...)

end # DrWatsonHelpers
