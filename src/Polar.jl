__precompile__()

module Polar

    using Requires

    include("Utils.jl")
    include("GF2n.jl")
    include("Stacks.jl")
    include("CommunicationsUtils.jl")

    include("PolarCodes.jl")
    include("PolarDecoding.jl")

    include("Distributions.jl")
    include("DensityEvolution.jl")
    include("Design.jl")

    include("Channels.jl")
    include("BridgeCh2DE.jl")

    include("SimulationUtils.jl")

    function __init__()
        @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("Visualizations.jl")
    end

end
