# 👑 Royal MD: the GOLD Standard in portable Molecular Dynamics
# Coded by Gleb Novikov
# ver 1.00: Royal Simulations with proteins, DNA, RNA and MORE ..
# last update 11/01/2026:
# - added support of all water models for AMBER/CHARMM force fields
# - added NVT/NPT equilibration with flexible restrains
# - added print_intro_banner to print intro logo with random pictures
# - improved output log in the production run (temperature, volume monitoring ..)
# (c) The VisualHub 2026

import os
import sys
import glob
import time
import random
from openmm.app import *
from openmm import *
from openmm.unit import *
from pdbfixer import PDBFixer


def print_intro_banner(add_slogan_separator=True):
    # Champagne Gold Palette
    GOLD_D  = "\033[38;5;172m" # Deep Gold (Bonds/Springs)
    GOLD_M  = "\033[38;5;220m" # Bright Gold (Atom Borders)
    GOLD_L  = "\033[38;5;221m" # Champagne Gold (Atom Centers)
    
    BOLD    = "\033[1m"
    RESET   = "\033[0m"

    # Here's the logo:
    logo = rf"""
{GOLD_M}{BOLD}
               Welcome to the Royal MD!

                 {GOLD_L}▟██▙{RESET}    {GOLD_L}▟██▙{RESET}    {GOLD_L}▟██▙{RESET}
                 {GOLD_L}▜██▛{RESET}    {GOLD_L}▜██▛{RESET}    {GOLD_L}▜██▛{RESET}
                  {GOLD_D}▜█▙{RESET}     {GOLD_D}█{RESET}     {GOLD_D}▟█▛{RESET}
                   {GOLD_M}▜██▙{GOLD_D}  {GOLD_M}▟█▙{GOLD_D}  {GOLD_M}▟██▛{RESET}
                    {GOLD_M}▜█▙{GOLD_D}▀ {GOLD_M}▜█▛{GOLD_D} ▀{GOLD_M}▜█▛{RESET}
                          {GOLD_D}█{RESET}
                          {GOLD_D}█{RESET}
                        {GOLD_L}▟██▙{RESET}
                        {GOLD_L}▜██▛{RESET}
{RESET}
    """
    print(logo)
    time.sleep(0.5)
    slogans = ["Automate", "Precise", "Beautiful"]
    separators = ["🥂", "👑", "🔱", "✨"]
    available_seps = separators.copy()
    print("          ", end="")
    # interate over slogans and print them in the banner
    for i, word in enumerate(slogans):
            time.sleep(1.0)
            
            if i == 0:
                # First word never needs a separator
                print(f"{BOLD}{GOLD_M}{word}", end="", flush=True)
            else:
                # Use specific boolean to decide on separator
                if add_slogan_separator and available_seps:
                    sep = random.choice(available_seps)
                    available_seps.remove(sep)
                    # create random brand logo
                    print(f" {sep} ", end="", flush=True)
                    print(word, end="", flush=True)
                else:
                    print(f"    {word}", end="", flush=True)
    time.sleep(2.0)
    print(RESET + "\n") # End the line after the loop

# --- STEP 0.1: CLEAN OLD DATA ---
def clean_old_files(target_pdb, activate: bool = False):
    if not activate:
        return
    print(f"🧹Cleaning Old Data:")
    patterns = ["*.dcd", "*.nc", target_pdb, "*.cif", "ligand_param/*"]
    for pattern in patterns:
        for file in glob.glob(pattern):
            try:
                if os.path.isfile(file):
                    os.remove(file)
                print(f"Removed: {file}")
            except Exception as e:
                pass

# --- STEP 1: FIX PROTEIN STRUCTURE ---
def fix_protein_topology(input_pdb, env_pH, ligands_to_remove):
    print(f"--- Step 1: Cleaning Protein Topology ---")
    fixer = PDBFixer(filename=input_pdb)

    # 1. STRIP LIGANDS FIRST (Stops the NoneType error)
    modeller = Modeller(fixer.topology, fixer.positions)
    residues_to_strip = [r for r in modeller.topology.residues() if r.name in ligands_to_remove]
    
    if residues_to_strip:
        removed_names = {r.name for r in residues_to_strip}
        print(f"🧹 Pre-cleaning: Removing {len(residues_to_strip)} residues {removed_names}")
        modeller.delete(residues_to_strip)
        
        # Synchronize fixer so it never "sees" the ligand
        fixer.topology = modeller.topology
        fixer.positions = modeller.positions
        # 👑 THE FIX: Clear the "memory" of the old gaps
        fixer.missingResidues = {}
    else:
        print("✨ No blacklisted residues found to strip.")

    # 2. Convert MSE to MET
    mse_count = 0
    cys_count = 0
    for residue in fixer.topology.residues():
        if residue.name == 'MSE':
            residue.name = 'MET'
            mse_count += 1
            for atom in residue.atoms():
                if atom.element.symbol == 'Se':
                    atom.element = Element.getBySymbol('S')
                    atom.name = 'SD' 
        elif residue.name in ['CYM', 'CYX']:
            residue.name = 'CYS'
            cys_count += 1


    if mse_count > 0:
        print(f"🧬 Converted {mse_count} 'MSE' residues to 'MET'.")
    if cys_count > 0:
        print(f"🧪 Standardized {cys_count} 'CYM/CYX' residues to 'CYS'.")

    # 3. Analyze and Repair (Only on clean protein)
    fixer.findMissingResidues()
    # added 29/12: quick GAP diagnostics (currently shows only gaps in the last chain ;-))
    if fixer.missingResidues:
        # Use [ ] to create a list explicitly, satisfying OpenMM's sum() 
        total_missing = sum([len(residues) for residues in fixer.missingResidues.values()])
        print(f"💫 GAPS DETECTED: Total of {total_missing} missing residues found.")
        
        chains = list(fixer.topology.chains())
        num_chains = len(chains)
        
        # Sort keys to ensure Chain A comes before Chain B
        sorted_keys = sorted(fixer.missingResidues.keys())
        
        for key in sorted_keys:
            chain_idx, res_idx = key
            residues = fixer.missingResidues[key]
            
            if chain_idx < num_chains:
                chain_id = chains[chain_idx].id
            else:
                chain_id = f"Unknown({chain_idx})"
            
            gap_size = len(residues)
            real_resid = res_idx +1
            print(f"  - Chain {chain_id}: Gap at residue {real_resid} ({gap_size} residues: {', '.join(residues)})")
    else:
        print("📜 NO GAPS DETECTED! The STRUCTURE is OKAY 👌")
    
    unique_residues_kept = {r.name for r in fixer.topology.residues()}
    print(f"⛩️ Residues kept: {sorted(list(unique_residues_kept))}")

    print(f"⚛️ Filling missing atoms and hydrogens at pH '{env_pH}'...")
    fixer.findMissingAtoms()
    fixer.addMissingAtoms() # Will not crash now
    fixer.addMissingHydrogens(env_pH)

    return fixer

# --- STEP 2: MERGE AND SOLVATE ---
def solvate_system(protein_fixer, ff_protein, ff_water, ff_map, box_type, box_padding, ionic_strength, pre_sim_pdb, input_pdb):
    print(f"-- Step 2: Solvation 🫧 --")

    # 1. Start with the fixed protein
    modeller = Modeller(protein_fixer.topology, protein_fixer.positions)

    # 2. Load Forcefields
    protein_xml = ff_map[ff_protein]
    try:
        water_xml = f"{ff_protein}/{ff_water}.xml"
        # added: test the ForceField creation the valided path
        _checkFF = ForceField(protein_xml, water_xml)
        ff_args = [protein_xml, water_xml]
    except:
        water_xml = f"{ff_water}.xml"
        ff_args = [protein_xml, water_xml]
    
    print(f"🥂 Successfully loading: {ff_args}")
    forcefield = ForceField(*ff_args)

    # --- FIX FOR 4-POINT WATER (like OPC with amber19) ---
    # .. check if we need the 4-point 'tip4pew' generator
    solvent_model = ff_water
    if 'opc' in ff_water.lower() or 'tip4p' in ff_water.lower():
        print(f"💧 4-point water detected ({ff_water}). Setting geometry to 'tip4pew'.")
        solvent_model = 'tip4pew'
    elif 'tip3p' in ff_water.lower() or ff_water.lower() == 'water':
        # CHARMM 'water' or standard 'tip3p' both use tip3p geometry
        print(f"💧 3-point water detected ({ff_water}). Setting geometry to 'tip3p'.")
        solvent_model = 'tip3p'
    elif 'spce' in ff_water.lower():
        print(f"💧 SPC/E water detected ({ff_water}). Setting geometry to 'spce'.")
        solvent_model = 'spce'
    else:
        print(f"😅 Unknown water '{ff_water}', defaulting to 'tip3p' geometry.")
        solvent_model = 'tip3p'

    # 3. Solvate
    print("💦 Adding solvent...")
    modeller.addSolvent(forcefield, 
                        padding=box_padding*nanometer, 
                        model=solvent_model,
                        boxShape=box_type, 
                        ionicStrength=ionic_strength*molar)
    
    print(f"Total System Size: {modeller.topology.getNumAtoms()} atoms.")
    with open(pre_sim_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
        
    return modeller, forcefield

# --- STEPS 3 - 5: MINIMIZE, EQUILIBRATE AND RUN ---
def setup_and_run(modeller, forcefield, temperature, pressure, timestep, equil_time, production_time, report_interval, output_nc, log_freq, pulling_config=None):
    print("--- Step 3: Setup & Minimize ⚡ ---")

    # Force conversion to numbers
    timestep = float(timestep)
    report_interval = int(report_interval)
    # production run steps:
    total_steps = int(production_time / timestep)

    # Create the System, Integrator, and Simulation
    system = forcefield.createSystem(modeller.topology, 
                                    nonbondedMethod=PME, 
                                    nonbondedCutoff=1.0*nanometer, 
                                    constraints=HBonds,
                                    rigidWater=True) # for OPC / TIP3P

    barostat = MonteCarloBarostat(pressure*bar, temperature*kelvin)
    system.addForce(barostat)

    # 3. Multi-step NPT equilibration with Harmonic Restraints
    equil_time = int(equil_time)
    num_stages = 10
    time_per_stage = equil_time / num_stages # for print report only
    k_max = 1000.0

    #formula = "0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)" # this may work only for vacum simulations
    formula = "0.5 * k * periodicdistance(x, y, z, x0, y0, z0)^2" # this works with PBC
    restraint_force = CustomExternalForce(formula)
    restraint_force.addGlobalParameter("k", k_max) 
    restraint_force.addPerParticleParameter("x0")
    restraint_force.addPerParticleParameter("y0")
    restraint_force.addPerParticleParameter("z0")

    protein_backbone = ['CA', 'C', 'N']
    nucleic_backbone = ['P', "O3'", "O5'", "C3'", "C4'", "C5'"]
    prot_count = 0
    na_count = 0

    # Get positions from modeller since simulation context doesn't exist yet
    positions_nm = modeller.positions.value_in_unit(nanometer)

    for atom in modeller.topology.atoms():
        res_name = atom.residue.name.strip()

        is_prot = atom.name in protein_backbone
        is_na = atom.name in nucleic_backbone

        if is_prot or is_na:
            pos = positions_nm[atom.index]
            restraint_force.addParticle(atom.index, [pos[0], pos[1], pos[2]])

            if is_prot: prot_count += 1
            if is_na: na_count += 1

    system.addForce(restraint_force)

    pull_head_force = None
    pull_tail_force = None
    pull_force_value = None

    def _residue_number(residue):
        try:
            return int(residue.id)
        except Exception:
            return residue.index + 1

    def _select_atoms_by_ranges(topology, ranges):
        indices = []
        for atom in topology.atoms():
            chain_id = atom.residue.chain.id
            resnum = _residue_number(atom.residue)
            for sel_chain, sel_start, sel_end in ranges:
                if chain_id == sel_chain and sel_start <= resnum <= sel_end:
                    indices.append(atom.index)
                    break
        return indices

    def _compute_com(atom_indices, positions_nm, system_obj):
        total_mass = 0.0
        com = Vec3(0.0, 0.0, 0.0)
        for idx in atom_indices:
            mass = system_obj.getParticleMass(idx).value_in_unit(dalton)
            if mass <= 0.0:
                continue
            pos = positions_nm[idx]
            com += Vec3(pos[0], pos[1], pos[2]) * mass
            total_mass += mass
        if total_mass == 0.0:
            raise ValueError("Selected atoms have zero total mass.")
        return com * (1.0 / total_mass)

    if pulling_config and pulling_config.get("enabled"):
        head_ranges = pulling_config.get("head_ranges", [])
        tail_ranges = pulling_config.get("tail_ranges", [])
        if not head_ranges or not tail_ranges:
            raise ValueError("Pulling is enabled but head/tail ranges are empty.")

        head_atoms = _select_atoms_by_ranges(modeller.topology, head_ranges)
        tail_atoms = _select_atoms_by_ranges(modeller.topology, tail_ranges)
        if not head_atoms or not tail_atoms:
            raise ValueError("Pulling selection produced empty head/tail atom groups.")

        head_com = _compute_com(head_atoms, positions_nm, system)
        tail_com = _compute_com(tail_atoms, positions_nm, system)
        direction = head_com - tail_com
        norm = (direction[0]**2 + direction[1]**2 + direction[2]**2) ** 0.5
        if norm == 0.0:
            raise ValueError("Head/tail COMs are identical; cannot define pull direction.")

        nx, ny, nz = direction[0] / norm, direction[1] / norm, direction[2] / norm
        pull_force_value = (pulling_config.get("force_pn", 20.0) * piconewton).value_in_unit(
            kilojoule_per_mole / nanometer
        )

        pull_head_force = CustomExternalForce("-f*(nx*x + ny*y + nz*z)")
        pull_head_force.addGlobalParameter("f", 0.0)
        pull_head_force.addGlobalParameter("nx", nx)
        pull_head_force.addGlobalParameter("ny", ny)
        pull_head_force.addGlobalParameter("nz", nz)
        for idx in head_atoms:
            pull_head_force.addParticle(idx, [])

        pull_tail_force = CustomExternalForce("f*(nx*x + ny*y + nz*z)")
        pull_tail_force.addGlobalParameter("f", 0.0)
        pull_tail_force.addGlobalParameter("nx", nx)
        pull_tail_force.addGlobalParameter("ny", ny)
        pull_tail_force.addGlobalParameter("nz", nz)
        for idx in tail_atoms:
            pull_tail_force.addParticle(idx, [])

        system.addForce(pull_head_force)
        system.addForce(pull_tail_force)
        print(f"🧲 Pulling enabled: head atoms={len(head_atoms)} tail atoms={len(tail_atoms)}")
    # Print the verification report
    print(f"🔒 Restraint Report:")
    print(f"   - Protein backbone: {prot_count}")
    print(f"   - DNA/RNA backbone: {na_count}")

    # note: changed collisions from 1 to 2/ps
    integrator = LangevinMiddleIntegrator(temperature*kelvin, 2/picosecond, timestep*picoseconds)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    # Detect Platform
    platform = simulation.context.getPlatform()
    print(f"✨ Running on Platform: {platform.getName()}")
    if platform.getName() in ['CUDA', 'OpenCL', 'Metal']:
        device_name = platform.getPropertyValue(simulation.context, 'DeviceName')
        print(f"💠 Active GPU: {device_name}")
    else:
        print("💻 Running on CPU (no GPU detected)")


    # --- STEP 3: THREE-STEP MINIMIZATION ---
    print("✨ Step 1/3: Minimizing solvent (Restraints HEAVY)... ", end="", flush=True)
    simulation.context.setParameter("k", k_max) 
    simulation.minimizeEnergy(tolerance=20.0*kilojoule_per_mole/nanometer)
    print("DONE!")

    print("✨ Step 2/3: Minimizing entire system (Restraints LIGHT)... ", end="", flush=True)
    # Don't go to 0.0 yet! 1.0 keeps the fold while letting sidechains pack.
    simulation.context.setParameter("k", 10.0)
    simulation.minimizeEnergy(tolerance=10.0*kilojoule_per_mole/nanometer)
    print("DONE!")

    print("✨ Step 3/3: Minimizing entire system (NO Restraints)... ", end="", flush=True)
    # Don't go to 0.0 yet! 1.0 keeps the fold while letting sidechains pack.
    simulation.context.setParameter("k", 0.0)
    simulation.minimizeEnergy(tolerance=2.0*kilojoule_per_mole/nanometer)
    print("DONE!")

    # --- DUMP MINIMIZED SNAPSHOT ---
    min_state = simulation.context.getState(getPositions=True, enforcePeriodicBox=False)
    current_positions = min_state.getPositions()
    #print(f"DEBUG: Topology atoms: {simulation.topology.getNumAtoms()} | State atoms: {len(current_positions)}")
    with open('minimized.cif', 'w') as f:
        PDBxFile.writeFile(simulation.topology, current_positions, f)

    # --- STEP 4: EQUILIBRATION (Chil Start) ---
    print(f"--- Step 4: System Equilibration 🌀 ---")

    simulation.context.setVelocitiesToTemperature(temperature*kelvin)

    # 4.1 NVT heating
    nvt_time = 100
    nvt_steps = int(nvt_time / timestep) # 100 ps

    # Disable Barostat
    barostat.setFrequency(0)
    simulation.context.reinitialize(preserveState=True)
    # Set High Restraints on backbone/ligands
    simulation.context.setParameter("k", k_max)

    # Run NVT equilibration to stabilize temperature first
    try:
        from mdtraj.reporters import NetCDFReporter
        simulation.reporters.append(NetCDFReporter('heating.nc', report_interval))
    except:
        simulation.reporters.append(DCDReporter('heating.dcd', report_interval))
    print(f"🌡️  Step 1/2: NVT equilibration ({nvt_time} ps):")
    print(f" 🔸 Heating system: T = {temperature:.1f}K ... ", end="", flush=True)
    simulation.step(nvt_steps)
    print("DONE!")

    # 4.2 Step-by-step NPT equilibration
    print(f"🧭 Step 2/2: NPT equilibration ({num_stages} × {time_per_stage:.0f} ps):")
    # Clear NVT reporter
    simulation.reporters.clear()
    # Turn on Barostat (default - every 25 steps)
    barostat.setFrequency(25)
    simulation.context.reinitialize(preserveState=True)

    try:
        from mdtraj.reporters import NetCDFReporter
        simulation.reporters.append(NetCDFReporter('equilibration.nc', report_interval))
    except:
        simulation.reporters.append(DCDReporter('equilibration.dcd', report_interval))

    # 3. GRADUATED EQUILIBRATION LOOP
    k_values = [k_max - (i * (k_max / (num_stages - 1))) for i in range(num_stages)]
    steps_per_stage = int((equil_time / timestep) / len(k_values))

    for k_val in k_values:
        print(f" 🔹 Relaxing system: k = {k_val:.1f}... ", end="", flush=True)
        simulation.context.setParameter("k", k_val)
        simulation.step(steps_per_stage)
        print("DONE!")

    print("👏 Equilibration finished!")

    # --- DUMP EQUILIBRATED SNAPSHOT ---
    eq_state = simulation.context.getState(getPositions=True, enforcePeriodicBox=False)
    eq_positions = eq_state.getPositions()
    with open('equilibrated.cif', 'w') as f:
        PDBxFile.writeFile(simulation.topology, eq_positions, f)

    # Turn off restraints for Production by setting force constant k to 0
    simulation.context.setParameter("k", 0.0)
    if pull_head_force and pull_tail_force:
        simulation.context.setParameter("f", pull_force_value)
    
    # Clear the equilibration reporter to start a fresh production file
    #simulation.reporters.pop()
    simulation.reporters.clear()

    # --- STEP 5: PRODUCTION ---
    # Add Production Reporters
    try:
        from mdtraj.reporters import NetCDFReporter
        simulation.reporters.append(NetCDFReporter(output_nc, report_interval))
    except:
        simulation.reporters.append(DCDReporter('production.dcd', report_interval))

    print(f"--- Step 5: Production Run 👑 ---")
    steps_per_percent = int(total_steps / 100)

    print(f"\n🔮 Simulation Duration: {total_steps * timestep / 1000:.2f} ns")
    print(f"\n{'%':>4} {'Step':>10} {'Pot. Energy':>15} {'Temp(K)':>8} {'Volume(nm³)':>8} {'Speed(ns/day)':>10}")
    print("-" * 55)

    production_start = time.time()
    last_time = time.time()

    for i in range(1, 101):
            simulation.step(steps_per_percent)
            
            if i % log_freq == 0:
                # FIX:
                state = simulation.context.getState(getEnergy=True)
                pot_energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
                kin_energy = state.getKineticEnergy()

                # GET VOLUME (To verify NPT)
                box_vectors = state.getPeriodicBoxVectors()
                volume = state.getPeriodicBoxVolume().value_in_unit(nanometer**3)
                
                # get the Temperature:
                # FIX for OPC water: Ccunt only atoms that have mass > 0
                real_atoms_dof = 0
                for atom in modeller.topology.atoms():
                    # Get mass from the system using the atom index
                    if system.getParticleMass(atom.index).value_in_unit(dalton) > 0:
                        real_atoms_dof += 3
                dof = real_atoms_dof - system.getNumConstraints()
                # Using OpenMM units for Boltzmann and Avogadro constants
                temp = (2 * kin_energy / (dof * BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA)).value_in_unit(kelvin)
                
                current_time = time.time()
                elapsed = current_time - last_time
                speed = (steps_per_percent * log_freq * timestep * 0.001) / (elapsed / 86400)
                
                print(f"{i:>3}% {simulation.currentStep:>10} {pot_energy:>15.1f} {temp:>8.1f} {volume:>10.1f} {speed:>10.1f}")
                last_time = current_time

    # Final Summary
    production_end = time.time()
    total_seconds = production_end - production_start
    avg_speed = (total_steps * timestep * 0.001) / (total_seconds / 86400)

    print("-" * 55)
    print(f"WORK COMPLETED!")
    print(f"⏳Total wall-clock time: {int(total_seconds // 60)}m {int(total_seconds % 60)}s")
    print(f"⚜️ Average performance: {avg_speed:.2f} ns/day")
    print(f"🌀Trajectory saved to: {output_nc}")

# --- MAIN CONTROLLER ---
def main():
    if len(sys.argv) < 2:
        print("ERROR: No input PDB file provided.")
        print("Usage: python royalMD.py perfetto.pdb")
        sys.exit(1)

    # --- MAIN CONFIG ---
    input_pdb = sys.argv[1]
    pre_sim_pdb = 'solvated.pdb'
    output_nc = 'production.nc'
    timestep = 0.002
    equil_time = 200 # 200ps of NPT equilibration (after 100ps of heating)
    production_time = 1000 # 1 ns of production
    report_interval = 1500 # report info every 1500 steps
    log_freq = 10
    # conditions
    temperature = 310 
    pressure = 1.0 
    # box and solvation: dodecahedron or octahedron
    box_type = 'dodecahedron'
    box_padding = 1.0 # if temperature fluctuates (with very flexible systems), increase it to 1.5
    ionic_strength = 0.15
    env_pH = 7.0
    # main params
    ff_protein = 'amber14' # NB: use amber14 with tip3p water
    ff_water = 'tip3p' # tip3p - with amber14/ amber99; opc - with amber19; water - with charmm36
    ligands_to_remove = ['SO4','HOH', 'EDO', 'LIH', 'LIG'] # .. add more from your PDBs
    # pulling config (constant force along initial head->tail vector)
    pulling_config = {
        "enabled": False,
        "force_pn": 20.0,  # 10-40 pN suggested
        # integrin avb3 example ranges (adjust to your topology)
        "head_ranges": [("A", 1, 500), ("B", 1, 500)],
        "tail_ranges": [("A", 900, 962), ("B", 650, 692)],
    }
    # NB: can be expanded with more force fields:
    ff_map = {
        'amber19': 'amber19-all.xml', # with OPC water
        'amber14': 'amber14-all.xml', # with TIP3P water
        'amber99sb': 'amber99sb.xml', # with TIP3P water
        'amber99sbildn': 'amber99sbildn.xml', # with TIP3P water
        'amber03': 'amber03.xml', # with TIP3P water
        'charmm36': 'charmm36.xml' # with "water" water ;-)
    }
    # --- 👑 EXECUTION PIPELINE 👑 ---
    print_intro_banner()
    # 0.1: Clean old trajectories, coordinate filles
    clean_old_files(pre_sim_pdb, activate=True)
    
    # 1. Fix Protein (stipping out ligands)
    fixer = fix_protein_topology(input_pdb, env_pH, ligands_to_remove)
    
    # 2: Solvate (Merging clean ligand inside)
    modeller, forcefield = solvate_system(fixer, ff_protein, ff_water, ff_map, 
                                          box_type, box_padding, ionic_strength, pre_sim_pdb, input_pdb)
    
    # 3 & 4: Run calculations
    setup_and_run(modeller, 
                  forcefield, 
                  temperature, 
                  pressure, 
                  timestep, 
                  equil_time,
                  production_time, 
                  report_interval, 
                  output_nc, 
                  log_freq,
                  pulling_config)

if __name__ == "__main__":
    main()
