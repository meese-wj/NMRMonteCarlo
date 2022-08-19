module DrWatsonHelpers

import DrWatson: scriptsdir

export slurmscriptsdir

slurmscriptsdir(args...) = scriptsdir("Slurm", args...)

end # DrWatsonHelpers