using OrbitalElements


const bc, M, G = 1.,1. ,1.
const numiso = NumericalIsochrone(G=G,M=M,bc=bc)
const anaiso = AnalyticalIsochrone(G=G,M=M,bc=bc)
const Ω0 = Ω₀(numiso)
const params = OrbitalParameters(Ω₀=Ω0)

println("Type of numiso:",typeof(numiso))

######################################################################
#
# Variables values
#
######################################################################

n1, n2 = -1, 2
a, e = 1.0, 0.5
u = 0.4

println("ANALYTICAL VALUES:")

E = OE.EFromAE(anaiso,a,e)
println("Energy: ",E)

J, L = OE.ComputeActionsAE(anaiso,a,e)
println("Actions: ",J," ; ",L)

aback, eback = OE.AEFromEL(anaiso,E,L)
println("Back from E, L: ",aback," ; ",eback)

Ω1, Ω2 = OE.ComputeFrequenciesAE(anaiso,a,e)
println("Frequencies: ",Ω1," ; ",Ω2)

aback, eback = OE.ComputeAEFromFrequencies(anaiso,Ω1,Ω2)
println("Back from frequencies: ",aback," ; ",eback)

println("NUMERICAL VALUES:")

E = OE.EFromAE(numiso,a,e,params)
println("Energy: ",E)

J, L = OE.ComputeActionsAE(numiso,a,e,params)
println("Actions: ",J," ; ",L)

aback, eback = OE.AEFromELBrute(E,L,numiso,params)
println("Back from E, L: ",aback," ; ",eback)

Ω1, Ω2 = OE.ComputeFrequenciesAE(numiso,a,e,params)
println("Frequencies: ",Ω1," ; ",Ω2)

aback, eback = OE.ComputeAEFromFrequencies(numiso,Ω1,Ω2,params)
println("Back from frequencies: ",aback," ; ",eback)
