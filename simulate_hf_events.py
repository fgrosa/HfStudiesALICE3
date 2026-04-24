#!/usr/bin/env python3

import sys
import os
import argparse
import yaml

from acts.examples.reconstruction import (
    addSeeding,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedFilterConfigArg,
    SpacePointGridConfigArg,
    SeedingAlgorithmConfigArg,
    SeedingAlgorithm,
    addCKFTracks,
    TrackSelectorConfig,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    CkfConfig,
    addVertexFitting,
    VertexFinder
)
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    ParticleConfig,
    addPythia8,
    addGenParticleSelection,
    addSimParticleSelection,
    addGeant4,
    addFatras,
    ParticleSelectorConfig,
    addDigitization,
)
import pathlib
import acts
import acts.examples
import acts.examples.geant4
import alice3

UNITS = acts.UnitConstants

def parse_pythia_cfg(config):
    """
    Docstring for parse_pythia_cfg
    
    --------------------------------
    PARAMETERS
    config: PYTHIA config file (.cfg)
    """

    system = acts.PdgParticle.eProton
    energy = 13600.
    args = []
    with open(config, "r") as f:
        for line in f:
            line_strip = line.strip()
            if "#" in line_strip:
                continue
            if "Beams:id" in line_strip:
                pdg_beam = int(line_strip.split("=")[-1])
                if pdg_beam == 2212:
                    system = acts.PdgParticle.eProton
                elif pdg_beam == 1000822080:
                    system = acts.PdgParticle.eLead
            elif "Beams:eCM" in line_strip:
                energy = float(line_strip.split("=")[-1])/1000
            elif line_strip != "":
                args.append(line_strip)

    return system, energy, args


def run_simulation(config_file):
    """
    Main function to run simulation

    --------------------------------
    PARAMETERS
    config_file: YAML config file
    """

    # read config file
    with open(config_file, "r") as file:
        cfg = yaml.safe_load(file)

    if not os.path.isdir(cfg["simulation"]["outputdir"]):
        os.makedirs(cfg["simulation"]["outputdir"])
    if not os.path.isdir(cfg["reconstruction"]["outputdir"]):
        os.makedirs(cfg["reconstruction"]["outputdir"])

    # we first build ALICE3 detector
    detector = alice3.buildALICE3Geometry(
        pathlib.Path("geom"),
        cfg["simulation"]["enable_material"],
        False,
        acts.logging.INFO
    )
    tracking_geom = detector.trackingGeometry()
    detector.contextDecorators()

    # initialise B-field
    if cfg["simulation"]["b_field_map"] != "":
        print("INFO: Using field map")
        field = acts.examples.MagneticFieldMapXyz(cfg["simulation"]["b_field_map"])
    else:
        print("INFO: Using constant B-field")
        field = acts.ConstantBField(acts.Vector3(0.0, 0.0, cfg["simulation"]["b_field"] * UNITS.T))

    # run generation
    cfg_sim = [cfg["simulation"]["pythia"]["enable"], cfg["simulation"]["gun"]["enable"], cfg["simulation"]["reader"]["enable"]]
    if sum(cfg_sim) > 1:
        print("ERROR: PYTHIA and gun generators cannot be enabled simultaneously, check your config file. Exit")
        sys.exit()
    if sum(cfg_sim) == 0:
        print("ERROR: One between PYTHIA and gun generators must be enabled, check your config file. Exit")
        sys.exit()

    rnd = acts.examples.RandomNumbers(seed=cfg["simulation"]["rnd_seed"])
    sequencer = acts.examples.Sequencer(events=cfg["simulation"]["n_events"],
                                        numThreads=cfg["simulation"]["n_threads"])

    if cfg["simulation"]["pythia"]["enable"]:
        pythia_cfg = parse_pythia_cfg(cfg["simulation"]["pythia"]["config"])
        pythia_cfg_pu = parse_pythia_cfg(cfg["simulation"]["pythia"]["config_pileup"])

        sequencer = addPythia8(
            sequencer,
            npileup=cfg["simulation"]["pythia"]["pileup"],
            beam=pythia_cfg[0],
            cmsEnergy=pythia_cfg[1] * acts.UnitConstants.TeV,
            hardProcess=pythia_cfg[2],
            pileupProcess=pythia_cfg_pu[2],
            vtxGen=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(
                    cfg["simulation"]["diamond"][0] * UNITS.mm,
                    cfg["simulation"]["diamond"][1] * UNITS.mm,
                    cfg["simulation"]["diamond"][2] * UNITS.mm,
                    0.0001 * UNITS.ns, # why?
                ),
                mean=acts.Vector4(0, 0, 0, 0),
            ),
            # printPythiaEventListing="long",
            rnd=rnd,
            logLevel=acts.logging.INFO,
            outputDirRoot=cfg["simulation"]["outputdir"],
            writeHepMC3=pathlib.Path(os.path.join(cfg["simulation"]["outputdir"], "events.hepmc3")),
            searchUpToHfQuark=cfg["simulation"]["pythia"]["search_hf_orig_up_to_quark"]
        )
    if cfg["simulation"]["gun"]["enable"]:
        # list of available particles: https://github.com/acts-project/acts/blob/main/Core/include/Acts/Definitions/PdgParticle.hpp#L21
        addParticleGun(
            sequencer,
            MomentumConfig(
                cfg["simulation"]["gun"]["pt"][0] * UNITS.GeV,
                cfg["simulation"]["gun"]["pt"][1] * UNITS.GeV,
                transverse=True
            ),
            EtaConfig(
                cfg["simulation"]["gun"]["eta"][0],
                cfg["simulation"]["gun"]["eta"][1],
                uniform=True
            ),
            ParticleConfig(
                cfg["simulation"]["gun"]["multiplicity"],
                acts.PdgParticle(cfg["simulation"]["gun"]["pdg"]),
                randomizeCharge=True,
            ),
            vtxGen=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(
                    cfg["simulation"]["diamond"][0] * UNITS.mm,
                    cfg["simulation"]["diamond"][1] * UNITS.mm,
                    cfg["simulation"]["diamond"][2] * UNITS.mm,
                    0.0001 * UNITS.ns, # why?
                ),
                mean=acts.Vector4(0, 0, 0, 0),
            ),
            rnd=rnd,
            logLevel=acts.logging.INFO,
            outputDirRoot=cfg["simulation"]["outputdir"],
        )
    elif cfg["simulation"]["reader"]["enable"]:
        sequencer.addReader(
            acts.examples.RootParticleReader(
                level=acts.logging.INFO,
                filePath=os.path.join(cfg["simulation"]["reader"]["inputdir"], "particles.root"),
                outputParticles="particles_generated",  # _generated",
                # outputSimHits="simhits",
                # particleCollection="particles",
                # inputDir="output",
                # inputFile="pythia8_particles.root",
            )
        )
        sequencer.addWhiteboardAlias("particles", "particles_generated")

    addGenParticleSelection(
        sequencer,
        ParticleSelectorConfig(
            eta=(-4., 4.),
            pt=(0.001 * UNITS.MeV, None),
            removeNeutral=False,
            # rho=(0.0, 24 * u.mm),
            # absZ=(0.0, 1.0 * u.m),
        ),
    )

    # perform transport
    if cfg["simulation"]["transport"] == "Geant4":
        gdml_detector = acts.examples.geant4.GdmlDetector(path="geom/o2sim_geometry.gdml")
        addGeant4(
            sequencer,
            gdml_detector,
            tracking_geom,
            field,
            materialMappings=["TRK_SILICON", "TF3_SILICON"],
            outputDirRoot=cfg["simulation"]["outputdir"],
            rnd=rnd,
            logLevel=acts.logging.DEBUG,  # acts.logging.INFO
            killVolume=tracking_geom.highestTrackingVolume,
            killAfterTime=40*UNITS.ns,  # 25 * u.ns,
        )
    elif cfg["simulation"]["transport"] == "Fatras":
        addFatras(
            sequencer,
            tracking_geom,
            field,
            enableInteractions=True,
            rnd=rnd,
            pMin=0.01,  # GeV, May 2025
            outputDirRoot=cfg["simulation"]["outputdir"],
            logLevel=acts.logging.INFO,
        )
    else:
        print(f"ERROR: transport option {cfg['simulation']['transport']} not supported. Exit")
        sys.exit()

    addSimParticleSelection(
        sequencer,
        ParticleSelectorConfig(
            hits=(7, None),
        ),
    )

    # digitisation
    sequencer = addDigitization(
        sequencer,
        tracking_geom,
        field,
        digiConfigFile=cfg["simulation"]["digi_file"],
        outputDirRoot=cfg["simulation"]["outputdir"],
        rnd=rnd,
        logLevel=acts.logging.DEBUG,
        # doMerge=True,
    )
    # for details on segmented digi: class Channelizer {..} https://github.com/acts-project/acts/blob/v36.3.2/Fatras/include/ActsFatras/Digitization/Channelizer.hpp#L21
    # auto channelsRes = m_channelizer.channelize() - https://github.com/acts-project/acts/blob/v36.2.1/Examples/Algorithms/Digitization/src/DigitizationAlgorithm.cpp#L212

    # reconstruction
    seeding_layers_opt = cfg["reconstruction"]["seeding"]["layers"]
    cfg_seeding_layers = None
    if seeding_layers_opt == "VD":
        cfg_seeding_layers = ("geom/geoSelectionForSeedingInner_BOTH_Barrel_Endcaps_NEW_VOLUMES_Feb2025.json")
    elif seeding_layers_opt == "ML3":
        cfg_seeding_layers = ("geom/geoSelectionForSeeding_3firstML.json")
    elif seeding_layers_opt == "MLall":
        cfg_seeding_layers = ("geom/geoSelectionForSeeding_5ML.json")
    elif seeding_layers_opt == "VDML":
        cfg_seeding_layers = ("geom/geoSelectionForSeeding_IB_ML.json")
    else:
        print(f"ERROR: seeding layers option {seeding_layers_opt} not supported. Exit")
        sys.exit()

    seeding_alg = SeedingAlgorithm.TruthSmeared
    if cfg["reconstruction"]["seeding"]["algorithm"] not in ["TruthSmeared", "GridTriplet"]:
        print(f"ERROR: seeding algorithm {cfg['reconstruction']['seeding']['algorithm']} not supported. Exit")
        sys.exit()
    elif cfg["reconstruction"]["seeding"]["algorithm"] == "GridTriplet":
        seeding_alg = SeedingAlgorithm.GridTriplet

    collision_region_4seeds = (cfg["reconstruction"]["seeding"]["collision_region"]) # mm; large values - for V0 daughter reconstruction

    # seeding
    sequencer = addSeeding(
        sequencer,
        tracking_geom,
        field,
        SeedFinderConfigArg(
            r=(None, 210 * UNITS.mm),  # iTOF is at 190 mm! if we want it for seeding
            deltaR=(1 * UNITS.mm, 200 * UNITS.mm),  # deltaR=(1. * u.mm, 60 * u.mm),
            collisionRegion=(-collision_region_4seeds * UNITS.mm, collision_region_4seeds * UNITS.mm),
            z=(-1000 * UNITS.mm, 1000 * UNITS.mm),
            maxSeedsPerSpM=cfg["reconstruction"]["seeding"]["max_seeds_per_spm"],  # 2 is minimum, >2 is better for Pb-Pb
            sigmaScattering=cfg["reconstruction"]["seeding"]["sigma_scattering"],
            # more info: https://github.com/acts-project/acts/blob/main/Core/include/Acts/Seeding/SeedFinderConfig.hpp
            radLengthPerSeed=cfg["reconstruction"]["seeding"]["radlen"],
            minPt=cfg["reconstruction"]["seeding"]["pt_min"] * UNITS.GeV,
            impactMax=cfg["reconstruction"]["seeding"]["imppar_max"] * UNITS.mm,  # important! IB vs ML seeds (e.g. 1 mm is ok for IB seeds, 5 mm - for ML seeds)
            cotThetaMax=cfg["reconstruction"]["seeding"]["costheta_max"],
            seedConfirmation=True,
            centralSeedConfirmationRange=acts.SeedConfirmationRangeConfig(
                zMinSeedConf=-620 * UNITS.mm,
                zMaxSeedConf=620 * UNITS.mm,
                rMaxSeedConf=4.9 * UNITS.mm,  # 36 * u.mm,  # IA: dramatically affects acceptance at eta ~4. <5 * u.mm  gives best results
                nTopForLargeR=1,  # number of top space points that confirm my seed at larger R, 1 - no confirmation
                nTopForSmallR=2,
            ),
            forwardSeedConfirmationRange=acts.SeedConfirmationRangeConfig(
                zMinSeedConf=-1220 * UNITS.mm,
                zMaxSeedConf=1220 * UNITS.mm,
                rMaxSeedConf=26 * UNITS.mm,  # 15 * u.mm,  #36 * u.mm,
                nTopForLargeR=1,
                nTopForSmallR=2,
            ),
            # skipPreviousTopSP=True,
            useVariableMiddleSPRange=True,
            deltaRMiddleSPRange=(0.2 * UNITS.mm, 1.0 * UNITS.mm),
        ),
        SeedFinderOptionsArg(bFieldInZ=cfg["simulation"]["b_field"] * UNITS.T, beamPos=(0 * UNITS.mm, 0 * UNITS.mm)),
        SeedFilterConfigArg(
            seedConfirmation=True if cfg["reconstruction"]["seeding"]["imppar_max"] < 2.0 else False,  # mm
            # If seedConfirmation is true we classify seeds as "high-quality" seeds.
            # Seeds that are not confirmed as "high-quality" are only selected if no
            # other "high-quality" seed has been found for that inner-middle doublet
            # Maximum number of normal seeds (not classified as "high-quality" seeds)
            # in seed confirmation
            maxSeedsPerSpMConf=1,  # 1 - USED_FOR_AUG_2025,#3, # CRUCIAL!!!!!!
            # 1 - USED_FOR_AUG_2025,   # Core/include/Acts/Seeding/SeedFilterConfig.hpp
            maxQualitySeedsPerSpMConf=1,
            # Maximum number of "high-quality" seeds for each inner-middle SP-dublet in
            # seed confirmation. If the limit is reached we check if there is a normal
            # quality seed to be replaced
        ),
        SpacePointGridConfigArg(
            impactMax=1.0 * UNITS.mm,
            phiBinDeflectionCoverage=3,
        ),
        SeedingAlgorithmConfigArg(
        ),
        # particleHypothesis=acts.ParticleHypothesis.pion, # IA
        # particleHypothesis=acts.ParticleHypothesis.electron, # IA
        geoSelectionConfigFile=cfg_seeding_layers,
        # seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
        seedingAlgorithm=seeding_alg,
        outputDirRoot=cfg["reconstruction"]["outputdir"],
        #    initialVarInflation = (50,50,50,50,50,50)  # IA
        #    initialVarInflation = (0.2,0.2,0.2,0.2,0.2,0.2)  #IA
    )

    #tracking
    # useful info: https://github.com/acts-project/acts/blob/main/Core/include/Acts/TrackFinding/MeasurementSelector.hpp
    sequencer = addCKFTracks(
        sequencer,
        tracking_geom,
        field,
        TrackSelectorConfig(
            pt=(0.06 * UNITS.GeV, 120 * UNITS.GeV),
            nMeasurementsMin=cfg["reconstruction"]["tracking"]["n_meas_min"],
            maxSharedHits=cfg["reconstruction"]["tracking"]["n_hits_shared_max"],
        ),  # IA, was: 500.0 * u.MeV
        ckfConfig=CkfConfig(
            chi2CutOffOutlier=cfg["reconstruction"]["tracking"]["chi2_outlier_max"],
            chi2CutOffMeasurement=cfg["reconstruction"]["tracking"]["chi2_meas_max"],
            numMeasurementsCutOff=cfg["reconstruction"]["tracking"]["meas_per_surf_max"],
            seedDeduplication=cfg["reconstruction"]["tracking"]["seed_deduplication"],
            stayOnSeed=cfg["reconstruction"]["tracking"]["stay_on_seed"],
        ),
        twoWay=cfg["reconstruction"]["tracking"]["two_way_ckf"],  # default: True,
        outputDirRoot=cfg["reconstruction"]["outputdir"],
        writeTrackSummary=False,
        logLevel=acts.logging.INFO,
    )

    sequencer = addAmbiguityResolution(
        sequencer,
        AmbiguityResolutionConfig(
            maximumSharedHits=cfg["reconstruction"]["ambiguity_resolution"]["n_hits_shared_max"],
            nMeasurementsMin=cfg["reconstruction"]["ambiguity_resolution"]["n_meas_min"]
        ),
        outputDirRoot=cfg["reconstruction"]["outputdir"],
        logLevel=acts.logging.INFO,
    )

    # vertexing
    sequencer = addVertexFitting(
        sequencer,
        field,
        vertexFinder=VertexFinder.AMVF,
        outputDirRoot=cfg["reconstruction"]["outputdir"],
        seeder=acts.examples.VertexSeedFinder.AdaptiveGridSeeder,
        useTime=False,  # True,
    )

    sequencer.run()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--config", "-c", metavar="text",
                        default="config.yml", help="YAML config file")
    args = parser.parse_args()

    run_simulation(args.config)
