# Simulate a single temperature
using DrWatson
@quickactivate :NMRMonteCarlo
using MonteCarloMeasurementUncertainty
import OnlineLogBinning: BinningAnalysisResult

sim = CleanNMRATMSimulation(; Lx = 8, βvalue = 0.2, Ntherm =2^18)
timer = @timed simulate!(sim)
@info "Simulation time: $(round(timer.time; sigdigits=4)) seconds"
mc_params = SimulationParameters(sim)
@info "$(round( (thermalization_sweeps(mc_params) + sampling_sweeps(mc_params)) * num_DoF(Hamiltonian(sim)) / timer.time; sigdigits = 4) ) updates/second"
@info "$( timer.bytes / 1000 ) KiB allocated"