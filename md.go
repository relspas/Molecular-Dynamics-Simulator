package main

import (
	"fmt"
)

//sigma and epsilon are parameters used by Lennard-Jones
type Particle struct {
	Position [3]float64
	Velocity [3]float64
	Mass     float64
	Sigma    float64
	Epsilon  float64
	Charge   float64
}

//constraint between two particles connected by a bond
type BondConstraint struct {
	Particle1Index int     // index of the first particle in the bond
	Particle2Index int     // index of the second particle in the bond
	SpringConstant float64 // Spring constant for the bond
	EquilibriumLen float64 // Equilibrium bond length
}

var ENABLE_OSCILLATION bool = false
var ENABLE_COULOMB bool = false

//Lennard-Jones potential, Coulomb, and harmonic behavior
func calculateForces(particles []Particle, constraints []BondConstraint) [][3]float64 {
	numParticles := len(particles)
	forces := make([][3]float64, numParticles)


	//Lennard-Jones force
	for i := 0; i < numParticles; i++ {
		for j := i + 1; j < numParticles; j++ {
			//Calculate the distance between particles i and j
			r_ij := distanceBetweenParticles(particles[i].Position, particles[j].Position)

			//Calculate Lennard-Jones potential and force for 2 particles
			epsilon_ij := math.Sqrt(particles[i].Epsilon * particles[j].Epsilon)
			sigma_ij := 0.5 * (particles[i].Sigma + particles[j].Sigma)

			//ljPotential := 4 * epsilon_ij * ((math.Pow(sigma_ij/r_ij, 12)) - (math.Pow(sigma_ij/r_ij, 6)))

			//the Lennard-Jones force between two particles at distance r
			// V(r) = 4 * epsilon * [(sigma/r)^12 - (sigma/r)^6]
			//source said "the force between two particles at distance r is the negative gradient of the potential is
			// F = -dV(r)/dr = (24 sigma^6 * epsilon * (r^6 - 2 * sigma^6))/r^13 (wolfram alpha derivation)
			ljForce := 24 * math.Pow(sigma_ij, 6) * epsilon_ij * (math.Pow(r_ij, 6) - 2*math.Pow(sigma_ij, 6)) / math.Pow(r_ij, 13)

			//direction of the force
			forceDirection := [3]float64{
				(particles[j].Position[0] - particles[i].Position[0]) / r_ij,
				(particles[j].Position[1] - particles[i].Position[1]) / r_ij,
				(particles[j].Position[2] - particles[i].Position[2]) / r_ij,
			}

			//Apply forces to particles i and j
			for k := 0; k < 3; k++ {
				forces[i][k] += ljForce * forceDirection[k]
				forces[j][k] -= ljForce * forceDirection[k]
			}
		}
	}


	// Coulomb forces between charged particles
	if ENABLE_COULOMB{
		for i := 0; i < numParticles; i++ {
			for j := i + 1; j < numParticles; j++ {
				qi := particles[i].Charge
				qj := particles[j].Charge

				// distance between particles i and j
				r_ij := distanceBetweenParticles(particles[i].Position, particles[j].Position)

				// Calculate Coulomb force
				// TODO need to add constant in but also need to figure out units
				// F = k*q1*q2/r^2
				coulombForce := (qi * qj) / (r_ij * r_ij) // Coulomb's law

				// Calculate the direction of the force
				forceDirection := [3]float64{
					(particles[j].Position[0] - particles[i].Position[0]) / r_ij,
					(particles[j].Position[1] - particles[i].Position[1]) / r_ij,
					(particles[j].Position[2] - particles[i].Position[2]) / r_ij,
				}

				// Apply Coulomb forces to particles i and j
				for k := 0; k < 3; k++ {
					forces[i][k] += coulombForce * forceDirection[k]
					forces[j][k] -= coulombForce * forceDirection[k]
				}
			}
		}
	}

	//Bond forces from harmonic behavior
	if ENABLE_OSCILLATION{
		for _, constraint := range constraints {
			p1 := particles[constraint.Particle1Index]
			p2 := particles[constraint.Particle2Index]
			bondVector := [3]float64{
				p2.Position[0] - p1.Position[0],
				p2.Position[1] - p1.Position[1],
				p2.Position[2] - p1.Position[2],
			}
			currentBondLength := math.Sqrt(bondVector[0]*bondVector[0] + bondVector[1]*bondVector[1] + bondVector[2]*bondVector[2])
			delta := currentBondLength - constraint.EquilibriumLen

			//hooke's Law (F = -k * delta)
			bondForce := constraint.SpringConstant * delta

			//Apply bond force in the direction of the bond vector to both particles
			for j := 0; j < 3; j++ {
				forces[constraint.Particle1Index][j] += bondForce * bondVector[j] / currentBondLength
				forces[constraint.Particle2Index][j] -= bondForce * bondVector[j] / currentBondLength
			}
		}
	}

	return forces
}

//SHAKE algorithm to enforce bond constraints
func applyConstraints(particles []Particle, constraints []BondConstraint, dt float64) {
	tolerance := 1e-5
	maxIterations := 100

	for iteration := 0; iteration < maxIterations; iteration++ {
		//Calculate the current bond lengths and differences from the desired lengths
		for _, constraint := range constraints {
			p1 := particles[constraint.Particle1Index]
			p2 := particles[constraint.Particle2Index]
			bondVector := [3]float64{
				p2.Position[0] - p1.Position[0],
				p2.Position[1] - p1.Position[1],
				p2.Position[2] - p1.Position[2],
			}
			currentBondLength := math.Sqrt(bondVector[0]*bondVector[0] + bondVector[1]*bondVector[1] + bondVector[2]*bondVector[2])
			delta := (currentBondLength - constraint.EquilibriumLen) / currentBondLength

			//Adjust particle positions to satisfy constraints
			for j := 0; j < 3; j++ {
				p1.Position[j] += 0.5 * delta * bondVector[j]
				p2.Position[j] -= 0.5 * delta * bondVector[j]
			}
		}

		//Recalculate forces after position adjustments
		forces := calculateForces(particles, constraints)

		//Check if constraints are satisfied
		constraintSatisfied := true
		for _, constraint := range constraints {
			p1 := particles[constraint.Particle1Index]
			p2 := particles[constraint.Particle2Index]
			bondVector := [3]float64{
				p2.Position[0] - p1.Position[0],
				p2.Position[1] - p1.Position[1],
				p2.Position[2] - p1.Position[2],
			}
			currentBondLength := math.Sqrt(bondVector[0]*bondVector[0] + bondVector[1]*bondVector[1] + bondVector[2]*bondVector[2])
			if math.Abs(currentBondLength-constraint.EquilibriumLen)/constraint.EquilibriumLen > tolerance {
				constraintSatisfied = false
				break
			}
		}

		if constraintSatisfied {
			break
		}

		updateVelocities(particles, forces, dt)
		updatePositions(particles, forces, dt)
	}
}


func distanceBetweenParticles(pos1, pos2 [3]float64) float64 {
    dx := pos2[0] - pos1[0]
    dy := pos2[1] - pos1[1]
    dz := pos2[2] - pos1[2]
    return math.Sqrt(dx*dx + dy*dy + dz*dz)
}

func updatePositions(particles []Particle, forces [][3]float64, dt float64) {
	for i := range particles {
		for j := 0; j < 3; j++ {
			particles[i].Position[j] += particles[i].Velocity[j] * dt
		}
	}
}

func updateVelocities(particles []Particle, forces [][3]float64, dt float64) {
	for i := range particles {
		for j := 0; j < 3; j++ {
			particles[i].Velocity[j] += forces[i][j] / particles[i].Mass * dt
		}
	}
}

//save particle data
func saveData(particles []Particle, step int) {
	//just print for demo
	fmt.Println(step, particles)
}


func main() {
	dt := 0.001          //Time step size for integration
	numSteps := 1000     //Number of simulation steps

	//Initialize two particles representing nitrogen and oxygen atoms
	particles := []Particle{
		{
			Position: [3]float64{0.0, 0.0, 0.0}, //Nitrogen atom at the origin
			Velocity: [3]float64{0.0, 0.0, 0.0},
			Mass:     14.0, //units?
			Sigma:    3.0,
			Epsilon:  1.0,
			Charge:   -0.5, //don't use this
		},
		{
			Position: [3]float64{1.0, 0.0, 0.0}, //Oxygen atom at (1.0, 0.0, 0.0)
			Velocity: [3]float64{0.0, 0.0, 0.0},
			Mass:     16.0, //units?
			Sigma:    3.5,
			Epsilon:  1.5,
			Charge:   0.3, //not using this
		},
	}

	//Add constraint for the covalent bond between nitrogen and oxygen

	constraints := []BondConstraint{
		{
			Particle1Index: 0, // Index of the oxygen atom in the particles slice
			Particle2Index: 1, // Index of the nitrogen      "        "
			SpringConstant: 1000.0, // bond stiffness
			EquilibriumLen: 1.2, // Desired bond length between O and N
		},
	}

	//Main simulation loop
	for step := 0; step < numSteps; step++ {
		//Apply constraints to maintain bond lengths
		applyConstraints(particles, constraints, dt)

		//Calculate forces based on the current positions and particle properties
		forces := calculateForces(particles, constraints)

		//Integrate equations of motion to update positions and velocities
		updatePositions(particles, forces, dt)
		updateVelocities(particles, forces, dt)

		//Output data and perform analysis
		if step%100 == 0 { //Change 100 to whatever output frequency
			saveData(particles, step)
		}
	}

	fmt.Println("Simulation completed.")
}
